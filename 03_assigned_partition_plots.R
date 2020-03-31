source('00_packages_functions.R', echo=TRUE)
#load.pigz(file = "flt3_aml_bakerlab.RData")

partition_assignment_plot <-
  custom_cp_plot(
    cds = cds_aligned_all,
    cp = "partition",
    alpha = 0.2,
    cell_size = 1,
    overwrite_labels = TRUE,
    plot_title = "Assigned Partitions",
    outfile = "plots_out/partition_assignment_plot.pdf",
    group_label_size = 5,
    w = 4.5,
    h = 4
  )

responses<-unique(colData(cds_aligned_all)$response)

pa_responses <- lapply(
  X = seq_along(responses),
  FUN = custom_cp_plot,
  cds = cds_aligned_all,
  cp = "partition",
  var = "response",
  value_to_highlight = responses,
  alpha = 0.2,
  cell_size = 1,
  overwrite_labels = TRUE,
  plot_title = responses,
  outfile = c("plots_out/par_refractory.pdf","plots_out/par_PR.pdf","plots_out/par_CR.pdf"),
  group_label_size = 5,
  w = 4.5,
  h = 4
)


colData(cds_aligned_all)$timepoint<-factor(colData(cds_aligned_all)$timepoint, levels = c("pretreatment","C1D26_gilt","post_1_cycle_gilt", "post_2_cycles_gilt","relapse"))
colData(cds_aligned_all)$pt_response<-paste0("Pt ",colData(cds_aligned_all)$patient,":  ",colData(cds_aligned_all)$response)
pa_pt_response_timepoint<-custom_variable_plot(cds_aligned_all,var = "partition_assignment",palette_viridis = F) + 
  facet_grid(rows = vars(pt_response), cols = vars(timepoint))+
  theme(strip.background = element_blank())
save_plot(pa_pt_response_timepoint, filename = "plots_out/pa_pt_response_timepoint.pdf", base_height = 5, base_width = 7.5)

cds_pt1_post2<-cds_aligned_all[,colData(cds_aligned_all)$sample_id %in% c("Y5209","Y6584")]
colData(cds_pt1_post2)
pt1_post2_tissue<-custom_variable_plot(cds_pt1_post2, var = "partition_assignment", palette_viridis = F) +
  facet_wrap(facets = "tissue")+
  theme(strip.background = element_blank())
save_plot(pt1_post2_tissue, filename = "plots_out/pt1_post2_tissue.pdf", base_width = 5, base_height = 2.5)


save.image.pigz(file = "flt3_aml_bakerlab.RData", n.cores = 39)
