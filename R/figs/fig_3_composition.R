fig_3 <- plot_grid(
  volcano_plot_pre,
  umap_bmx_response_timepoint,
  NULL,
  NULL,
  nrow = 2,
  align = "h", 
  axis = "b",
  rel_heights = c(1, 2),
  rel_widths = c(0.8, 1),
  labels = c("A", "B")
  
)


save_plot(
  plot = fig_3,
  # filename = "test.png",
  filename = str_glue("{network_out}/fig_3.png"),
  base_width = 7.5, 
  base_height = 9.75
)
