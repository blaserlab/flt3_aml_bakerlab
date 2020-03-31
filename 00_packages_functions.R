#libraries
library("monocle3")
library("tidyverse")
library("cowplot")
theme_set(theme_cowplot(font_size = 11))
library("fastSave")
library("ggrepel")
library("htmltools")
library("ggpubr")
library("pheatmap")
library("ggplotify")

#custom operators
`%notin%` <- Negate(`%in%`)

# custom functions

plot_cells_alt <-
  function (cds,
            gene_or_genes,
            h,
            w,
            cell_size = 1,
            alpha = 1 ,
            outfile = NULL,
            ncol = NULL,
            plot_title = NULL,
            plot_type = "pdf",
            i =1) {
    
    if (i>1) {
      outfile<-paste0("plots_out/",names(marker_gene_list)[i],".",plot_type)
      gene_or_genes<-gene_or_genes[[i]]
      plot_title<-plot_title[[i]]
    }
    data <- plot_cells(cds = cds, genes = gene_or_genes)[["data"]]
    data$gene_short_name <-
      factor(data$gene_short_name, levels = gene_or_genes)
    background_data <- data %>% filter(is.na(value))
    foreground_data <- data %>% filter(!is.na(value))
    p <- ggplot() +
      geom_point(
        data = background_data,
        aes(x = data_dim_1, y = data_dim_2),
        color = "grey80",
        shape = 1,
        size = cell_size,
        stroke = 0.25
      ) +
      geom_point(
        data = foreground_data,
        aes(
          x = data_dim_1,
          y = data_dim_2,
          color = log10(value)
        ),
        shape = 16,
        size = cell_size,
        alpha = alpha
      ) +
      scale_color_viridis_c(end = 0.8) +
      labs(
        x = "UMAP 1",
        y = "UMAP 2",
        color = "Log10(Expr)",
        title = plot_title
      ) +
      facet_wrap( ~ gene_short_name, ncol = ncol) +
      theme(strip.background = element_blank())+
      theme(plot.title = element_text(hjust = 0.5))
    if (!is.null(outfile)) {
      save_plot(
        plot = p,
        filename = outfile,
        base_height = h,
        base_width = w
      )
    }
    return(p)
  }


cumulative_max_expr<-function(extracted_df,gene_list){
  return(extracted_df %>% 
           tbl_df() %>% 
           select(-feature_id,-gene_short_name,-id) %>% 
           pivot_wider(names_from = feature_label, values_from = value) %>% 
           replace(., is.na(.),0) %>%
           mutate(max_val = do.call(pmax, c(select(., one_of(gene_list))))) %>%
           mutate(max_val = na_if(max_val,0)))
}

add_cds_factor_columns<-function(cds, columns_to_add){
  for (i in 1:length(columns_to_add)) {
    colData(cds)$new<-unname(columns_to_add[i])
    names(colData(cds))[names(colData(cds)) == "new"] <- names(columns_to_add[i])
  }
  return(cds)
}




custom_variable_plot<-function(cds,
                               i = NULL, 
                               var, 
                               value_to_highlight = NULL, 
                               foreground_alpha = 1, 
                               legend_pos = "right", 
                               cell_size = 0.5, 
                               legend_title = NULL, 
                               plot_title = NULL, 
                               outfile = NULL, 
                               h = NULL, 
                               w = NULL, 
                               palette_viridis = T) {
  if (!is.null(i)) {
    value_to_highlight<-value_to_highlight[[i]]
    plot_title<-plot_title[[i]]
    outfile<-paste0("plots_out/",outfile[[i]],".pdf")
  }
  
  data<-plot_cells(cds)[["data"]]
  data_long<-data %>% pivot_longer(cols = matches(var), names_to = "var")
  plot<-ggplot()
  if(!is.null(value_to_highlight)){
    data_background<-data_long %>% filter(value %notin% value_to_highlight)
    data_long<-data_long %>% filter(value %in% value_to_highlight)
    plot<-plot+
      geom_point(data = data_background, 
                 aes(x = data_dim_1, y = data_dim_2), 
                 stroke = 0.25, 
                 shape = 1, 
                 size = cell_size, 
                 color = "grey80")
  }
  plot<-plot+
    geom_point(data = data_long, 
               aes(x = data_dim_1, y = data_dim_2, fill = value, color = value), 
               stroke = 0.25, shape = 21, 
               alpha = foreground_alpha, 
               size = cell_size)
  if(class(data_long$value)=="numeric") {
    plot<-plot+scale_fill_viridis_c(aesthetics = c("color", "fill"))
  } else if (palette_viridis == T) {
    plot<-plot+
      scale_fill_viridis_d(begin = 0.1,end = 0.9)+
      scale_color_viridis_d(begin = 0.1,end = 0.9, guide = F)
  } else {
    plot<-plot+scale_color_discrete(guide = F)+
      scale_fill_discrete()
  }
  plot<-plot+guides(fill = guide_legend(override.aes = list(size = 2, alpha = 1, color = "transparent")))
  plot<-plot+labs(x = "UMAP 1", y = "UMAP 2", title = plot_title, fill = legend_title) + theme(plot.title = element_text(hjust = 0.5))
  plot<-plot+theme(legend.position = legend_pos)#+coord_fixed()
  if(is.null(outfile)) return(plot) else save_plot(plot = plot, filename = outfile, base_height = h, base_width = w);return(plot)
  
}

#test<-custom_cp_plot(cds = cds_aligned, cp = "cluster", var = "specimen", value_to_highlight = "specimen_0");test

custom_cp_plot <- function(cds,
                           i = 1,
                           var = NULL,
                           alpha = 1,
                           cp = c("cluster", "partition"),
                           overwrite_labels = T,
                           legend_pos = "none",
                           cell_size = 0.5,
                           legend_title = NULL,
                           plot_title = NULL,
                           outfile = NULL,
                           h = NULL,
                           w = NULL,
                           group_label_size = 3,
                           value_to_highlight = NULL,
                           plot_type = "pdf") {
  value_to_highlight<-value_to_highlight[[i]]
  plot_title<-plot_title[[i]]
  outfile<-outfile[[i]]
  #extract the data from the input cds
  data <- plot_cells(cds)[["data"]]
  #convert to long format
  data_long <- pivot_longer(data = data,
                            cols = cp,
                            names_to = "cp")
  # generate text data frame for cluster/partition labels
  text_df <- data_long %>% group_by(value)
  median_coord_df <-
    data_long %>% group_by(value) %>% summarise(
      fraction_of_group = n(),
      text_x = median(data_dim_1),
      text_y = median(data_dim_2)
    )
  text_df <- left_join(text_df, median_coord_df)
  if (cp == "cluster")
    text_df$label <-
    text_df$cluster_assignment
  else
    text_df$label <- text_df$partition_assignment
  text_df <-
    text_df %>% group_by(label) %>% summarise(text_x = dplyr::first(text_x),
                                              text_y = dplyr::first(text_y))
  # if highlighting a categorical variable, generate background data and keep data_long as foreground
  if (!is.null(var)) {
    background_data_long <- data_long %>% filter((!!sym(var)) != value_to_highlight)
    data_long <- data_long %>% filter((!!sym(var)) == value_to_highlight)
  }
  #initialize the plot
  plot <- ggplot()
  # lay down the background plot if using
  if(!is.null(value_to_highlight)){
    plot<-plot+
      geom_point(data = background_data_long,
                 aes(x = data_dim_1, y = data_dim_2),
                 stroke = 0.25,
                 shape = 1,
                 size = cell_size,
                 color = "grey80")
  }
  
  # make the main colored plot from data_long
  plot <- plot +
   geom_point(
      data = data_long,
      aes(
        x = data_dim_1,
        y = data_dim_2,
        fill = value,
        color = value),
      stroke = 0.25,
      shape = 21,
      alpha = alpha,
      size = cell_size) +
    scale_color_discrete(guide = F) +
    labs(
      x = "UMAP 1",
      y = "UMAP 2",
      title = plot_title,
      fill = legend_title) +
    theme(plot.title = element_text(hjust = 0.5))
  # overwrite labels if you want to
  if (overwrite_labels == T) {
    plot <- plot +
      theme(legend.position = "none") +
      ggrepel::geom_text_repel(
        data = text_df,
        mapping = aes_string(x = "text_x", y = "text_y", label = "label"),
        size = group_label_size) +
      guides(fill = guide_legend(override.aes = list(
        size = 2,
        alpha = 1,
        color = "transparent"
      )))
  } else {
    plot <- plot +
      theme(legend.position = legend_pos) +
      guides(fill = guide_legend(override.aes = list(
        size = 2,
        alpha = 1,
        color = "transparent"
      )))}
  
  # option to save plot
  if (is.null(outfile)) {return(plot)} else {
    save_plot(
      plot = plot,
      filename = outfile,
      base_height = h,
      base_width = w)
  return(plot)}
}

custom_violin_plot <-
  function(cds,
           variable,
           genes_to_plot,
           outfile = NULL,
           pseudocount = 1,
           include_jitter = FALSE,
           ytitle = "Log10(Expr+1)",
           plot_title = NULL,
           w,
           h,
           rows = 1,
           show_x_label = TRUE,
           legend_pos = "none"#,
           #comparison_list = NULL,
           #sig_lab_y = 1,
           #yplotmax
           ) {
    my_comparisons <-
      comparison_list#(list(c(comparator1,comparator2),c(comparator1,comparator3)...))
    data_to_plot <-
      plot_genes_violin(cds_subset = cds[rowData(cds)$gene_short_name %in% genes_to_plot,], group_cells_by = variable)[["data"]]
    p1 <-
      ggplot(data = data_to_plot, aes(
        x = !!as.name(variable),
        y = log10(expression +
                    pseudocount)
      )) #expression already normalized when data extracted by violin plot function
    p1 <- p1 +
      geom_violin(
        scale = "width",
        color = "black",
        trim = T,
        size = 0.5,
        aes(fill = !!as.name(variable)),
        draw_quantiles = 0.5
      )
    if (include_jitter == TRUE) {
      p1 <-
        p1 + geom_jitter(
          shape = 16,
          size = 0.05,
          color = "black",
          alpha = 0.1,
          width = 0.2
        )
    }
    # p1 <- p1 +
    #   ylim(0,yplotmax)
    # if (!is.null(comparison_list)) {
    #   p1<-p1+stat_compare_means(
    #     comparisons = my_comparisons,
    #     method = "wilcox.test",
    #     size = 2,
    #     label = "p.signif",
    #     hide.ns = F,
    #     label.y = sig_lab_y
    #   )
    # } 
    p1<-p1+
      theme(legend.position = legend_pos) +
      theme(legend.direction = "horizontal") +
      theme(legend.justification = "center") +
      #scale_y_continuous(breaks = unname(labels_breaks_vec), labels = names(labels_breaks_vec))+
      scale_fill_viridis_d(alpha = 0.6,
                           begin = 0.1,
                           end = 0.9) +
      labs(
        x = "",
        y = ytitle,
        title = plot_title,
        fill = NULL
      ) +
      theme(plot.title = element_text(hjust = 0.5)) +
      #theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
      facet_wrap(~ gene_short_name, nrow = rows) +
      theme(strip.background = element_rect(fill = "transparent"))
    if (show_x_label == F) {
      p1 <- p1 + theme(axis.text.x = element_blank())
    }
    if (!is.null(outfile)) {
      save_plot(p1, filename = outfile, base_width = w, base_height = h)
    }
    return(p1)
  }
