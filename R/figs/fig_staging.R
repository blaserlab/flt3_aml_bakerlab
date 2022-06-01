# partition umap -----------------------------------
umap_partition <- 
  bb_var_umap(
    cds_anno_aligned_tissue_id,
    var = "partition_assignment_1",
    palette = experimental_group_palette,
    cell_size = 0.5,
    foreground_alpha = 0.1,
    overwrite_labels = T
    
  )

# local cells umap faceted ----------------------------- 
umap_myeloblast_localn <- 
  bb_var_umap(
  cds_anno_aligned_tissue_id[, colData(cds_anno_aligned_tissue_id)$partition_assignment_1 == "Myeloblast"],
  var = "log_local_n",
  facet_by = c("timepoint_binary", "chemo_response"),
  rows = vars(chemo_response),
  cols = vars(timepoint_binary),
  scales = "free",
  cell_size = 0.5,
  foreground_alpha = 0.4,
  nbin = 150 
) + 
  theme(plot.title = element_text(hjust = 0.5)) +
  theme(panel.background = element_rect(color = "grey80")) +
  scale_x_continuous(n.breaks = 4) +
  labs(color = "Log<sub>10</sub><br>Local<br>Cells") +
  theme(legend.title = element_markdown()) +
  ggtitle("Myeloblasts")



# volcano pre -------------------------------------------------------------
genes_to_highlight_pre <- c("BMX", "JAK3")

volcano_data_pre <-
  pseudobulk_pre_myeloblast %>%
  mutate(threshold = padj < 0.1 & abs(log2FoldChange) >= 0.58) %>%
  mutate(text_label = ifelse(gene_short_name %in% genes_to_highlight_pre, gene_short_name, ""))

volcano_plot_pre <-
  ggplot(
    volcano_data_pre,
    aes(
      x = log2FoldChange,
      y = -log10(padj),
      colour = threshold,
      fill = threshold,
      label = text_label
    )
  ) +
  geom_point(shape = 21,
             size = 0.5,
             alpha = 0.4) +
  geom_text_repel(color = "black",
                  fontface = "italic",
                  box.padding = 0.5,
                  point.padding = 0.25,
                  min.segment.length = 0,
                  max.overlaps = 20000,
                  size = 3,
                  segment.size = 0.25,
                  force = 2,
                  seed = 1234,
                  segment.curvature = -0.1,
                  segment.square = TRUE,
                  segment.inflect = TRUE) +
  xlab("log<sub>2</sub> fold change") +
  ylab("-log<sub>10</sub> adjusted p-value") +
  theme(axis.title.x =  element_markdown()) +
  theme(axis.title.y = element_markdown()) +
  theme(legend.position = "none") +
  scale_color_manual(values = c("grey80", "#DC0000")) +
  scale_fill_manual(values = c("transparent", "#DC0000")) +
  labs(caption = "\U21D0 Up in Sensitive\nUp in Unresponsive \U21D2", title = "Pre-Treatment")+
  theme(plot.caption.position = "panel") +
  theme(plot.caption = element_text(hjust = 0.5)) +
  theme(plot.title = element_text(hjust = 0.5)) +
  coord_cartesian(xlim = c(-1.1*max(abs(range(volcano_data_pre %>% dplyr::filter(!is.na(padj)) %>% pull(log2FoldChange)))), 1.1*max(abs(range(volcano_data_pre %>% filter(!is.na(padj)) %>% pull(log2FoldChange))))))


# volcano post  -----------------------------------------------------------

genes_to_highlight_post <- c("CCL5", "CXCL1", "CXCL2")

volcano_data_post <- pseudobulk_post_myeloblast %>% 
  mutate(threshold = padj < 0.1 & abs(log2FoldChange) >= 0.58) %>%
  mutate(text_label = ifelse(gene_short_name %in% genes_to_highlight_post, gene_short_name, NA))

volcano_plot_post <-
  ggplot(
    volcano_data_post,
    aes(
      x = log2FoldChange,
      y = -log10(padj),
      colour = threshold,
      fill = threshold,
      label = text_label
    )
  ) +
  geom_point(shape = 21, 
             size = 0.5, 
             alpha = 0.4) +
  geom_text_repel(color = "black",
                  fontface = "italic",  
                  box.padding = 0.5,
                  point.padding = 0.25,
                  min.segment.length = 0,
                  max.overlaps = 2000,
                  nudge_x = -0.5,
                  size = 3, 
                  segment.size = 0.25,
                  force = 2,
                  seed = 1234,
                  segment.curvature = -0.1,
                  segment.square = TRUE,
                  segment.inflect = TRUE) +
  xlab("log<sub>2</sub> fold change") +
  ylab("-log<sub>10</sub> adjusted p-value") +
  theme(axis.title.x =  element_markdown()) +
  theme(axis.title.y = element_markdown()) +
  theme(legend.position = "none") +
  scale_color_manual(values = c("grey80", "#DC0000")) +
  scale_fill_manual(values = c("transparent", "#DC0000")) +
  labs(caption = "\U21D0 Up in Sensitive\nUp in Unresponsive \U21D2", title = "Post-Treatment")+
  theme(plot.caption.position = "panel") +
  theme(plot.caption = element_text(hjust = 0.5)) +
  theme(plot.title = element_text(hjust = 0.5)) + 
  coord_cartesian(xlim = c(-1.1*max(abs(range(volcano_data_post %>% filter(!is.na(padj)) %>% pull(log2FoldChange)))), 1.1*max(abs(range(volcano_data_post %>% filter(!is.na(padj)) %>% pull(log2FoldChange))))))
volcano_plot_post

# gene dotplot --------------------------------------------------------
dotplot_genes <-
  c(
    "FLT3",
    "IL3RA",
    "IL6R",
    "CSF2RA",
    "CXCL8",
    "CXCL2",
    "CCL5",
    "CCL3",
    "theend"
  )

dotplot_genes_revision <- CellChat::CellChatDB.human$interaction |> 
  filter(ligand %in% dotplot_genes) |> 
  bind_rows(CellChat::CellChatDB.human$interaction |> 
              filter(receptor %in% dotplot_genes)) |> 
  select(ligand, receptor, pathway_name) |> 
  as_tibble() |> 
  bind_rows(tribble(~ligand, ~receptor, ~pathway_name,
                    "IL6", "IL6R", "IL6",
                    "CSF2", "CSF2RA", "CSF2")) |> 
  pivot_longer(c(ligand, receptor), values_to = "Gene", names_to = "gene_type") |> 
  distinct() |> 
  filter(Gene != "ACKR1")

colData(cds_anno_aligned_tissue_id)$timepoint_binary_short <- recode(colData(cds_anno_aligned_tissue_id)$timepoint_binary, 
                                                                     "pre-treatment" = "pre",
                                                                     "post-treatment" = "post")

myeloblast_dotplot <- 
  bb_gene_dotplot(
  cds = cds_anno_aligned_tissue_id[,colData(cds_anno_aligned_tissue_id)$partition_assignment_1 == "Myeloblast"],
  # markers = dotplot_genes,
  markers = unique(dotplot_genes_revision$Gene),
  scale_expression_by_gene = TRUE,
  group_cells_by = "multifactorial",
  group_ordering = tribble(
  ~aesthetic,~variable, ~value, ~level,
  "facet","chemo_response", "sensitive", 1,
  "facet","chemo_response", "unresponsive", 2,
  "axis", "timepoint_binary_short", "pre", 1,
  "axis", "timepoint_binary_short", "post", 2
  ),
colorscale_name = "Expression",
sizescale_name = "Proportion\nExpressing",
# gene_ordering = dotplot_genes, 
max.size = 5
  
) + labs(x = NULL, y = NULL) 

myeloblast_dotplot_revision_data <- left_join(myeloblast_dotplot[["data"]], dotplot_genes_revision) |> 
  mutate(Gene = factor(Gene)) |> 
  mutate(gene_type = factor(gene_type, levels = c("receptor", "ligand"))) |> 
  mutate(pathway_name = factor(pathway_name, levels = rev(c("FLT3", "IL6", "CSF2", "CXCL", "CCL")))) |> 
  mutate(Gene = fct_reorder2(Gene, pathway_name, gene_type)) |> 
  mutate(Gene = fct_rev(Gene))

myeloblast_dotplot_revision <- ggplot(myeloblast_dotplot_revision_data, 
       mapping = aes(x = axis, y = Gene, color = mean, size = percentage)) +
  geom_point() +
  scale_color_viridis_c() +
  scale_size_area(max_size = 4) +
  facet_grid(pathway_name ~ facet, scales = "free", space = "free_y") +
  theme(strip.background = element_blank()) +
  theme(strip.text.y = element_blank()) +
  theme(axis.text.y = element_text(size = 8, face = "italic"))+
  labs(color = "Expression", 
       size = "Proportion\nExpressing",
       x = NULL,
       y = NULL)

# gene module heatmap ------------------------------------------------------

# set up the annotation dataframe
module_ha_df <-
  tibble(columns = colnames(agg_mat_blasts)) %>%
  mutate(Response = recode(columns, "sensitive_post-treatment" = "sensitive",
                                   "sensitive_pre-treatment" = "sensitive",
                                   "unresponsive_post-treatment" = "unresponsive",
                                   "unresponsive_pre-treatment" = "unresponsive")) %>%
  mutate(Timepoint = recode(columns, "sensitive_post-treatment" = "post-treatment",
                                   "sensitive_pre-treatment" = "pre-treatment",
                                   "unresponsive_post-treatment" = "post-treatment",
                                   "unresponsive_pre-treatment" = "pre-treatment")) %>%
  mutate(Timepoint = fct_relevel(Timepoint, c("pre-treatment", "post-treatment"))) %>%
  as.data.frame()
rownames(module_ha_df) <- module_ha_df$columns
module_ha_df <- module_ha_df[,-1]


col_fun_2 <-
  colorRamp2(breaks = c(min(scale(agg_mat_blasts)),
                        0,
                        max(scale(agg_mat_blasts))),
             colors = c("#3C5488", "white", "#DC0000"))

module_ha <- HeatmapAnnotation(
  df = module_ha_df,
  col = list(
    Response = c("sensitive" = "#7FBC41", "unresponsive" = "#F1B6DA"),
    Timepoint = c("pre-treatment" = "#9E0142", "post-treatment" = "#FEE08B")
  ),
  gp = gpar(col = "grey80"),
  annotation_legend_param = list(
    Response = list(
      ncol = 1,
      at = c("unresponsive", "sensitive")
    ),
    Timepoint = list(
      ncol = 1
    )
  ),
  show_annotation_name = T,
  annotation_name_gp = gpar(fontsize = 9)
)

gene_module_heatmap_blasts <-
  grid.grabExpr(draw(
    Heatmap(
      matrix = scale(agg_mat_blasts),
      name = "Module\nExpression",
      col = col_fun_2,
      show_row_dend = F,
      show_column_names = F,
      show_row_names = T,
      row_names_gp = gpar(fontsize = 8),
      top_annotation = module_ha,
      border = "grey80",
      rect_gp = gpar(col = "grey80"),
      column_dend_height = unit(3, "mm"),
      heatmap_legend_param = list(
        direction = "horizontal", 
        title_position = "lefttop"
      )
    ),padding = unit(c(10, 10, 5, 20), "mm"), 
    heatmap_legend_side = "bottom", 
    annotation_legend_side = "bottom"),wrap = T# bottom left top right
  )




# go term bubbles----------------------------------------------------

revigo_bubbles <-
  function(revigo_data_tbl,
           dispensability_thresh,
           pval_thresh,
           plot_text_size = 3,
           box_padding = 0.5,
           bubble_size = c(2, 5)) {
    revigo_data <- revigo_data_tbl %>%
      filter(plot_X != "null") %>%
      filter(plot_Y != "null") %>%
      mutate(
        plot_X = as.numeric(plot_X),
        plot_Y = as.numeric(plot_Y),
        plot_size = as.numeric(plot_size),
        log10_p = as.numeric(`log10 p-value`),
        frequency_pct = as.numeric(str_replace(frequency, "%", "")),
        uniqueness = as.numeric(uniqueness),
        dispensability = as.numeric(dispensability),
        eliminiated = recode(eliminated, "0" = FALSE, "1" = TRUE)
      ) %>%
      select(-frequency, -`log10 p-value`)
    text_data <- revigo_data %>%
      filter(dispensability <= dispensability_thresh) %>%
      filter(log10_p < pval_thresh) %>%
      filter(eliminated == F)
    p <- ggplot(revigo_data) +
      geom_point(
        aes(
          x = plot_X,
          y = plot_Y,
          fill = -log10_p,
          size = plot_size
        ),
        alpha = I(0.4),
        shape = 21
      ) +
      scale_size(range = bubble_size, guide = F) +
      scale_fill_viridis_c() +
      geom_label_repel(
        data = text_data,
        aes(x = plot_X, y = plot_Y, label = description),
        color = "black",
        box.padding =  box_padding,
        point.padding = 0,
        min.segment.length = 0,
        max.overlaps = 200,
        nudge_x = -0.5,
        size = plot_text_size,
        segment.size = 0.25,
        force = 1,
        seed = 1234,
        segment.curvature = -0.1,
        segment.square = TRUE,
        segment.inflect = TRUE,
        fill = alpha(c("white"),0.5)
      ) +
      labs(
        x = "semantic space x",
        y = "semantic space y",
        size = "GO term frequency",
        fill = "-log10 p"
      ) +
      theme(legend.key = element_blank())
    
    return(p)
  }


# edit the revigo data for conciseness-------------------------

revigo_data_trimmed <-
  revigo_data %>%
  mutate(description = recode(description,
                              "SRP-dependent cotranslational protein targeting to membrane" = 
                                "SRP-dependent cotranslational\nprotein targeting to membrane",
                              "regulation of multicellular organismal process" = 
                                "regulation of multicellular\norganismal process",
                              "mRNA transport" = "mRNA\ntransport",
                              "macromolecular complex subunit organization" = "macromolecular complex",
                              "cellular component organization or biogenesis" = "cellular component"))



mod1_bubbles <- revigo_bubbles(
  revigo_data_tbl = revigo_data_trimmed %>% filter(module == "module_1"),
  dispensability_thresh = 0.2,
  pval_thresh = -3,
  plot_text_size = 3,
  box_padding = 0.75,
  bubble_size = c(1, 3)
) +
  labs(title = "Module 1") +
  theme(plot.title = element_text(hjust = 0.5)) 

mod2_bubbles <- revigo_bubbles(
  revigo_data_tbl = revigo_data_trimmed %>% filter(module == "module_2"),
  dispensability_thresh = 0.06,
  pval_thresh = -3,
  plot_text_size = 3,
  box_padding = 0.25,
  bubble_size = c(1, 3)
) +  
  labs(title = "Module 2") +
  theme(plot.title = element_text(hjust = 0.5)) 

mod3_bubbles <- revigo_bubbles(
  revigo_data_tbl = revigo_data_trimmed %>% filter(module == "module_3"),
  dispensability_thresh = 0.08,
  pval_thresh = -17,
  plot_text_size = 3,
  box_padding = 0.5,
  bubble_size = c(1,3)
) +  
  labs(title = "Module 3") +
  theme(plot.title = element_text(hjust = 0.5))

mod4_bubbles <- revigo_bubbles(
  revigo_data_tbl = revigo_data_trimmed %>% filter(module == "module_4"),
  dispensability_thresh = 0.1,
  pval_thresh = -20,
  plot_text_size = 3,
  box_padding = 0.5,
  bubble_size = c(1, 3)
) +  
  labs(title = "Module 4") +
  theme(plot.title = element_text(hjust = 0.5)) 

mod5_bubbles <- revigo_bubbles(
  revigo_data_tbl = revigo_data_trimmed %>% filter(module == "module_5"),
  dispensability_thresh = 0.1,
  pval_thresh = -3,
  plot_text_size = 3,
  box_padding = 0.5,
  bubble_size = c(2, 6)
) +  
  labs(title = "Module 5") +
  theme(plot.title = element_text(hjust = 0.5)) 

mod13_bubbles <- revigo_bubbles(
  revigo_data_tbl = revigo_data_trimmed %>% filter(module == "module_13"),
  dispensability_thresh = 0.1,
  pval_thresh = -50,
  plot_text_size = 3,
  box_padding = 0.5,
  bubble_size = c(1,3)
) +  
  labs(title = "Module 13") +
  theme(plot.title = element_text(hjust = 0.5)) 
mod13_bubbles


# cluster membership heatmap -------------------------------------------- 

cluster_memb_mat <- 
  colData(cds_anno_aligned_tissue_id) %>%
  as_tibble() %>%
  group_by(sample_id,partition_assignment_1) %>%
  summarise(n = n()) %>%
  ungroup(partition_assignment_1) %>%
  mutate(freq = n/sum(n)) %>%
  select(-n) %>%
  pivot_wider(names_from = partition_assignment_1, values_from = freq, values_fill = 0) %>%
  bb_tbl_to_matrix() %>%
  t()


heatmap_metadata <- 
  colData(cds_anno_aligned_tissue_id) %>%
  as_tibble() %>%
  group_by(sample_id, tissue, chemo_response, timepoint_binary,gilt_paper_id_corrected) %>%
  summarise(n = n()) %>%
  ungroup() %>%
  select(sample_id, 
         Patient = gilt_paper_id_corrected, 
         Tissue = tissue, 
         Response = chemo_response, 
         Timepoint = timepoint_binary, 
         Cells = n) %>%
  as.data.frame()
rownames(heatmap_metadata) <- heatmap_metadata$sample_id
heatmap_metadata <- heatmap_metadata[,-1]

# color functions
col_fun <- colorRamp2(breaks = c(0,1), colors = c("transparent","#DC0000"))
col_fun_1 <- colorRamp2(breaks = c(0,max(heatmap_metadata$Cells)), colors = c("transparent", "darkgreen"))

ha <- HeatmapAnnotation(
  df = heatmap_metadata, # has to be in the same order as the columns it is labeling; there is no other connection
  col = list(Patient = experimental_group_palette,
             Tissue = experimental_group_palette,
             Response = experimental_group_palette,
             Timepoint = experimental_group_palette,
             Cells = col_fun_1
  ),
  gp = gpar(col = "grey80"),
  annotation_name_gp = gpar(fontsize = 10)
)

cluster_memb_heatmap <-
  grid.grabExpr(draw(
    ComplexHeatmap::Heatmap(
      matrix = cluster_memb_mat,
      col = col_fun,
      rect_gp = gpar(color = "grey80"),
      top_annotation = ha,
      column_title = "Samples",
      column_title_side = "top",
      column_dend_side = "top",
      row_names_gp = gpar(fontsize = 10),
      show_column_names = F,
      column_dend_height = unit(5,"mm"),
      row_dend_width = unit(5,"mm"),
      name = "Sample\nProportion"
    ), 
    padding = unit(c(25, 2, 2, 10), "mm"),# bottom left top right
    annotation_legend_side = "bottom",adjust_annotation_extension = TRUE 
    ),wrap = TRUE
    )

# canonical gene dotplots ------------------------------

canonical_gene_dotplots <- bb_gene_dotplot(cds = cds_anno_aligned_tissue_id, 
                markers = c("CD34",
                            "CD3E",
                            "CD14",
                            "CD79A",
                            "FCGR3A",
                            "HLA-DRA",
                            "LILRA4",
                            "MZB1","KIT","CD33"),
                group_cells_by = "partition_assignment_1",
                colorscale_name = "Expression", 
                sizescale_name = "Proportion\nExpressing") + 
  labs(x = NULL, y = NULL) +
  theme(axis.text.x =  element_text(angle = 30, hjust = 1))

# bmx umap ---------------------------------------------------

umap_bmx_response_timepoint <-
  bb_gene_umap(
    cds_anno_aligned_tissue_id,
    gene_or_genes = "BMX",
    cell_size = 1,
    alpha = 0.8,
    color_legend_title = "Expression"
  ) +
  facet_grid(cols = vars(timepoint_binary),
             rows = vars(chemo_response)) +
  labs(title = "BMX") +
  theme(plot.title = element_text(hjust = 0.5, face = "italic")) +
  theme(panel.background = element_rect(color = "grey80"))

# cytokine data ----------------------------------------------

cytokine_plotfun <- function(type, samp, title, dat, pal, yt) {
  dat_avg <- dat |>
    filter(value_type == "average") |>
    filter(sample == samp) |>
    filter(sample_type == type) |>
    filter(analyte != "CXCL1") |>
    filter(!is.na(value))
  dat_rep <- dat |>
    filter(str_detect(value_type, "rep")) |>
    filter(sample == samp) |>
    filter(sample_type == type) |>
    filter(analyte != "CXCL1") |>
    filter(!is.na(value))
  stat.test <- dat_rep %>%
    group_by(analyte) |> 
    t_test(value ~ condition) |>  
    add_significance()
  stat.test <- stat.test |> 
    add_xy_position(x = "analyte", 
                    dodge = 0.5, 
                    y.trans = yt) |> 
    filter(group1 %in% c("WT", "DMSO"))
  if ("CCL4" %notin% stat.test$analyte) {
    stat.test$x <- ifelse(stat.test$analyte %in% c("CCL5", "CXCL8"), 
                          stat.test$x+1, 
                          stat.test$x)
    stat.test$xmin <- ifelse(stat.test$analyte %in% c("CCL5", "CXCL8"), 
                          stat.test$xmin+1, 
                          stat.test$xmin)
    stat.test$xmax <- ifelse(stat.test$analyte %in% c("CCL5", "CXCL8"), 
                          stat.test$xmax+1, 
                          stat.test$xmax)
  }
  conds <- as.character(unique(dat_avg$condition))
  pal <- pal[conds]
  p <- ggplot(mapping = aes(
    x = analyte,
    y = value,
  )) +
    geom_bar(
      data = dat_avg,
      stat = "identity",
      width = 0.4,
      alpha = 0.2,
      position = position_dodge(width = 0.5),
      mapping = aes(
        fill = condition,
        color = condition
      )
    ) +
    geom_point(
      data = dat_rep,
      pch = 21,
      size = 2,
      position = position_jitterdodge(jitter.width = 0.2, dodge.width = 0.5),
      show.legend = FALSE,
      mapping = aes(
        fill = condition,
        color = condition
      )
    ) +
    ggpubr::stat_pvalue_manual(
      stat.test,
      # label = "p.adj",
      tip.length = 0.01, 
      bracket.nudge.y = 0.1, 
      hide.ns = TRUE
    ) +
    scale_color_manual(values = pal) +
    scale_fill_manual(values = alpha(colour = pal,
                                     alpha = 0.4)) +
    theme(legend.position = c(0.1, 0.9)) +
    theme(legend.title = element_blank()) +
    labs(x = NULL, y = "pg/mL", title = title) +
    scale_x_discrete(limits = c("CCL2", "CCL3", "CCL4", "CCL5", "CXCL8"))
  if (samp == "MOLM13") {
    p <- p +
      scale_y_log10() +
      annotation_logticks(sides = "l")
    
  }
  p
}



cytokine_plotlist <- pmap(
  .l = list(
    x = c("MV411", "MV411", "MOLM13"),
    y = c("bmx_crispr", "cell_line_treatment", "cell_line_treatment"),
    z = c("", "MV4-11", "MOLM13"),
    ytrans = c(\(x){x}, \(x){x}, log10)
  ),
  .f = \(x, y, z, ytrans) cytokine_plotfun(
    samp = x,
    type = y,
    title = z,
    dat = cytokine_array_data,
    pal = experimental_group_palette,
    yt = ytrans
  )
) |> set_names(c("crispr_ko_plot", "mv411_drug_plot", "molm13_drug_plot"))

# figure 5D
dmso_matrix <- 
  cytokine_array_data |> 
  filter(sample_type == "primary") |> 
  filter(value_type == "average") |> 
  filter(condition == "DMSO") |> 
  select(sample, value, analyte) |> 
  pivot_wider(names_from = analyte, values_from = value) |> 
  bb_tbl_to_matrix()

chmfl_matrix <- 
  cytokine_array_data |> 
  filter(sample_type == "primary") |> 
  filter(value_type == "average") |> 
  filter(condition == "CHMFL") |> 
  select(sample, value, analyte) |> 
  pivot_wider(names_from = analyte, values_from = value) |> 
  bb_tbl_to_matrix()


plot_matrix <- t(log2(chmfl_matrix / dmso_matrix))

col_fun <- circlize::colorRamp2(breaks = c(-4, 0, 4),
                                colors = c("red3", "white", "blue4"))
pt_sample_heatmap <-
  grid.grabExpr(draw(Heatmap(
    plot_matrix, 
    col = col_fun, name = "Log2 Fold Change", 
    row_dend_width = unit(5, "mm"), 
    column_dend_height = unit(5, "mm"), 
    row_names_gp = gpar(fontsize = 9), 
    column_names_gp = gpar(fontsize = 9),
    heatmap_legend_param = list(direction = "horizontal", 
                                position = "lefttop")
  ), heatmap_legend_side = "bottom"), wrap = T)
plot_grid(pt_sample_heatmap)
