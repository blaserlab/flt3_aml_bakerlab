possibly_subsetCommunication <- possibly(CellChat::subsetCommunication, otherwise = NULL)
blast_autocrine_mat <- map2_dfr(.x = cellchat_list,
         .y = names(cellchat_list),
         .f = \(x, y) {
           dat <- possibly_subsetCommunication(x) 
           if (!is.null(dat)) {
             dat <- dat |> 
             filter(source == "Myeloblast") |> 
             filter(target == "Myeloblast") |> 
             filter(pval < 0.05) |> 
             select(interaction_name, prob) |> 
             mutate(sample = y) |> 
             relocate(sample)
           }
    }) |> 
  pivot_wider(names_from = sample, values_from = prob, values_fill = 0) |> 
  mutate(interaction_name = str_replace_all(interaction_name, "_", "-")) |> 
  bb_tbl_to_matrix()

blast_autocrine_mat_col_fun <- circlize::colorRamp2(breaks = c(0, max(blast_autocrine_mat)), 
                                                    colors = c("transparent", "red4"))

blast_autocrine_col_anno_df <- tibble(colnames = colnames(blast_autocrine_mat)) |> 
  separate(col = colnames, c("Response", "Timepoint", "Patient"), sep = "_", remove = FALSE) |> 
  mutate(Patient = str_remove(Patient, "pt")) |> 
  mutate(Patient = factor(Patient, levels = c("2", "4", "5", "17", "20", "21", "22"))) |>
  mutate(Timepoint = factor(Timepoint, levels = c("pre-treatment", "post-treatment"))) |> 
  tibble::column_to_rownames(var = "colnames")

blast_autocrine_col_anno_fun <-
  HeatmapAnnotation(
    df = blast_autocrine_col_anno_df,
    col = list(Response = experimental_group_palette, 
               Timepoint = experimental_group_palette,
               Patient = experimental_group_palette),
    gp = gpar(col = "white"),
    annotation_name_gp = gpar(fontsize = 10, fontface = "bold")
  )

cellchat_myeloblast_heatmap <-
  grid.grabExpr(draw(
    Heatmap(
      blast_autocrine_mat,
      name = "Interaction\nScore",
      col = blast_autocrine_mat_col_fun,
      top_annotation = blast_autocrine_col_anno_fun,
      show_column_names = FALSE,
      cluster_columns = FALSE,
      column_split = blast_autocrine_col_anno_df$Response,
      border = TRUE,
      row_names_gp = gpar(fontsize = 10),
      column_title = "Inferred Myeloblast Autocrine Interactions",
      row_title = "Interaction",
      row_dend_width = unit(3, "mm"), 
      column_title_gp = gpar(fontsize = 12), 
      row_title_gp = gpar(fontsize = 12)
    )
  ),
  wrap = TRUE)

