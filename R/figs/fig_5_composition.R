
fig5_left <- plot_grid(
  cytokine_plotlist$crispr_ko_plot,
  cytokine_plotlist$mv411_drug_plot,
  cytokine_plotlist$ molm13_drug_plot,
  nrow = 4,
  align = "v",
  axis = "l",
  labels = c("A", "B", "C"),
  rel_heights = c(1, 1, 1, 1)
)

fig5_right <- plot_grid(
  pt_sample_heatmap,
  NULL,
  nrow = 2,
  labels = c("D", ""),
  rel_heights = c(1,1)
)

fig5 <- plot_grid(
  fig5_left,
  fig5_right,
  ncol = 2,
  rel_widths = c(2,1)
)

save_plot(
  plot = fig5+ theme(plot.background = element_rect(fill = "white", color = "white")), 
  filename = str_glue("{network_out}/fig_5.png"),
  # filename = "test.png",
  base_width = 7.5, 
  base_height = 9.75
)

