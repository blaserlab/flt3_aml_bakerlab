fig_S2 <- 
  plot_grid(
    cluster_memb_heatmap,
    canonical_gene_dotplots,
    nrow = 2,
    rel_heights = c(1.4,1),
    labels = c("A","B"),label_y = c(1,1.05)
  )

save_plot(
  plot = fig_S2 + theme(plot.background = element_rect(fill = "white", color = "white")),
  # filename = "test.png",
  filename = str_glue("{network_out}/fig_S2.png"),
  base_width = 7.5, 
  base_height = 9.75
)