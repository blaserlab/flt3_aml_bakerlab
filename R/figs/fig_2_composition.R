
fig2 <- plot_grid(
  umap_partition,
  umap_myeloblast_localn,
  volcano_plot_post,
  myeloblast_dotplot,
  NULL,
  NULL,
  nrow = 3,
  align = "vh",
  axis = "lb",
  labels = c("A", "B", "C", "D", "", ""),
  rel_heights = c(1, 1, 1)
)

save_plot(
  plot = fig2,
  # filename = "test.png",
  filename = str_glue("{network_out}/fig_2.png"),
  base_width = 7.5, 
  base_height = 9.75
)

