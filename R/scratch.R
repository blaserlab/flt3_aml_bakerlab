bb_gene_dotplot(
  filter_cds(
    cds_anno_aligned_tissue_id,
    cells = bb_cellmeta(cds_anno_aligned_tissue_id) |> filter(partition_assignment %in% c("Blast 1"))
  ),
  markers = c("IRF3", "IRF7", "CXCL8"),
  group_cells_by = "multifactorial",
  group_ordering = tribble(
  ~aesthetic,~variable, ~value, ~level,
  "facet","chemo_response", "sensitive", 1,
  "facet","chemo_response", "unresponsive", 2,
  "axis", "timepoint_binary", "pre-treatment", 1,
  "axis", "timepoint_binary", "post-treatment", 2
  ),
  scale_expression_by_gene = TRUE
    
  
) +
  theme(axis.text.x = element_text(angle = 30,hjust = 1)) +
  labs(x = NULL, y = NULL)

bb_var_umap(cds_anno_aligned_tissue_id, "partition_assignment", overwrite_labels = T)
bb_cellmeta(filter_cds(
  cds_anno_aligned_tissue_id,
  cells = bb_cellmeta(cds_anno_aligned_tissue_id) |> filter(partition_assignment %in% c("Blast 1"))
)) |> group_by(chemo_response, timepoint_binary) |> 
  summarise()
bb_gene_umap(cds_anno_aligned_tissue_id, "CXCL8") + facet_grid(chemo_response ~ timepoint_binary)
bb_gene_umap(cds_anno_aligned_tissue_id, "IRF3") + facet_grid(chemo_response ~ timepoint_binary)
