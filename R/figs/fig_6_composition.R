
fig_6_top_left <- 
  plot_grid(
    gene_module_heatmap_blasts,
    labels = "A"
  )

fig_6_top_right <- 
  plot_grid(
    mod2_bubbles,
    mod4_bubbles,
    nrow = 2,
    labels = c("B", "C"),
    align = "v", 
    axis = "l"
  )

fig_6_top <- 
  plot_grid(
    fig_6_top_left,
    fig_6_top_right,
    ncol = 2
  )

fig_6_bottom <- 
  plot_grid(
    mod13_bubbles,
    mod1_bubbles,
    mod3_bubbles,
    NULL,
    nrow = 2,
    labels = c("D","E","F")
    
  )


fig_6 <- 
  plot_grid(
    fig_6_top,
    fig_6_bottom,
    nrow = 2,
    rel_heights = c(1,1)
  )

cowplot::save_plot(
  plot = fig_6,
  # filename = "test.pdf",
  filename = str_glue("{network_out}/fig_6.pdf"),
  base_width = 7.5, 
  base_height = 9.75,
) 
