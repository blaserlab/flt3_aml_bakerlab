
fig_2.5 <- 
  plot_grid(
    cellchat_myeloblast_heatmap,
    NULL,
    nrow = 2
  )

cowplot::save_plot(
  plot = fig_2.5 + theme(plot.background = element_rect(fill = "white", color = "white")),
  # filename = "test.png",
  filename = file.path(network_out, "fig_S2.5.png"),
  base_width = 7.5, 
  base_height = 9.75,
) 
