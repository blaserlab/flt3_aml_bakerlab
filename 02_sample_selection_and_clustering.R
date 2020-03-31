source('00_packages_functions.R', echo=TRUE)
#load.pigz(file = "flt3_aml_bakerlab.RData")

# choose whether to include blood and whether to use alignment or not
# cds_trimmed_all - not aligned!
colData(cds_trimmed_all)$patient <-
  factor(colData(cds_trimmed_all)$patient)# fix cds patient data type
cds_trimmed_all_cvp <-
  custom_variable_plot(
    cds = cds_trimmed_all,
    var = "tissue",
    palette_viridis = TRUE,
    foreground_alpha = 0.2,
    cell_size = 1,
    plot_title = "Not Aligned",
    legend_title = "Tissue",
    outfile = "plots_out/cds_trimmed_all_cvp.pdf",
    w = 5,
    h = 4
  )

cds_trimmed_all_cvp_faceted <- cds_trimmed_all_cvp +
  facet_wrap(facets = "patient") +
  theme(strip.background = element_blank())
save_plot(
  cds_trimmed_all_cvp_faceted,
  filename = "plots_out/cds_trimmed_all_cvp_faceted.pdf",
  base_width = 5,
  base_height = 4
)

#cds aligned by patient variable
colData(cds_aligned_all)$patient <-
  factor(colData(cds_aligned_all)$patient)# fix cds patient data type
cds_aligned_all_cvp <-
  custom_variable_plot(
    cds = cds_aligned_all,
    var = "tissue",
    palette_viridis = TRUE,
    foreground_alpha = 0.2,
    cell_size = 1,
    plot_title = "Aligned by Patient",
    legend_title = "Tissue",
    outfile = "plots_out/cds_aligned_all_cvp.pdf",
    w = 5,
    h = 4
  )

cds_aligned_all_cvp_faceted <- cds_aligned_all_cvp +
  facet_wrap(facets = "patient") +
  theme(strip.background = element_blank())
save_plot(
  cds_aligned_all_cvp_faceted,
  filename = "plots_out/cds_aligned_all_cvp_faceted.pdf",
  base_width = 5,
  base_height = 4
)

# Group cells into clusters
cds_aligned_all <-
  cluster_cells(cds_aligned_all, cluster_method = "louvain")

colData(cds_aligned_all)$cluster <- clusters(cds_aligned_all)
colData(cds_aligned_all)$cluster_assignment <-
  clusters(cds_aligned_all)

colData(cds_aligned_all)$partition <- partitions(cds_aligned_all)
colData(cds_aligned_all)$partition_assignment <-
  partitions(cds_aligned_all)

#plot clusters and partitions prior to assignment

cluster_plot <-
  custom_cp_plot(
    cds = cds_aligned_all,
    cp = "cluster",
    alpha = 0.2,
    cell_size = 1,
    overwrite_labels = TRUE,
    plot_title = "Louvain Clusters",
    outfile = "plots_out/cluster_plot.pdf",
    group_label_size = 5,
    w = 4.5,
    h = 4
  )

partition_plot <-
  custom_cp_plot(
    cds = cds_aligned_all,
    cp = "partition",
    alpha = 0.2,
    cell_size = 1,
    overwrite_labels = TRUE,
    plot_title = "Partitions",
    outfile = "plots_out/partition_plot.pdf",
    group_label_size = 5,
    w = 4.5,
    h = 4
  )

# identify marker genes
marker_test_res_c <-
  top_markers(
    cds_aligned_all,
    group_cells_by = "cluster",
    reference_cells = 1000,
    cores = 39,
    genes_to_test_per_group = 50
  )

marker_test_res_p <-
  top_markers(
    cds_aligned_all,
    group_cells_by = "partition",
    reference_cells = 1000,
    cores = 39,
    genes_to_test_per_group = 50
  )

top_specific_markers <- marker_test_res_p %>%
  filter(fraction_expressing >= 0.10) %>%
  group_by(cell_group) %>%
  top_n(3, pseudo_R2)

top_specific_marker_ids <- unique(top_specific_markers %>% pull(gene_short_name))
top_specific_marker_ids<-top_specific_marker_ids[str_sub(top_specific_marker_ids,1,2)!="MT"]


#manually assign celltypes to clusters and partitions
cds_aligned_all$partition_assignment <- recode(
  cds_aligned_all$partition,
  "1" = "B1",#MS4A1 PAX5, IGHM
  "2" = "AML/Progenitor", #CSF3R, MSI2, GNA15
  "3" = "T1",#CD3E,CD6,CD247
  "4" = "PC",#XBP1, IGKC, MZB1
  "5" = "Erythroid",#GYPA, KLF1, GATA1 
  "6" = "T/NK",#IL32, TRAC, PRF1
  "7" = "MK/Plt",#CD41, GP9, PF4
  "8" = "Blast 1", #CD34,CD38,FLT3
  "9" = "Blast 2", #CXCR4,IL3RA,NAMPT
  "10" = "Mono/Mac", #CD68, LYZ, MPEG1
  "11" = "Blast 3", #MPO, ELANE, NPM1
)
#export partition top markers
left_join(
  marker_test_res_p,
  tbl_df(colData(cds_aligned_all)) %>% select(partition, partition_assignment) %>% distinct(),
  by = c("cell_group" = "partition")
) %>% arrange(partition_assignment) %>%
  write_csv("data_out/partition_markers.csv")

#export cluster top markers
marker_test_res_c %>% rename("cluster" = "cell_group") %>% write_csv("data_out/cluster_markers.csv")

# make marker gene plots
marker_gene_list<-list(
  B1_genes=c("MS4A1", "PAX5", "IGHM"),
  AML_prog_genes=c("CSF3R", "MSI2", "GNA15"),
  T1_genes=c("CD3E","CD6","CD247"),
  PC_genes=c("XBP1", "IGKC", "MZB1"),
  ery_genes=c("GYPA", "KLF1", "GATA1"),
  TNK_genes=c("IL32", "TRAC", "PRF1"),
  MK_plt_genes=c("CD41", "GP9", "PF4"),
  bl_1_genes=c("CD34","CD38","FLT3"),
  bl_2_genes=c("CXCR4","IL3RA","NAMPT"),
  monomac_genes=c("CD68", "LYZ", "MPEG1"),
  bl_3_genes=c("MPO", "ELANE", "NPM1")
  
)

titles<-tbl_df(colData(cds_aligned_all)) %>% 
  pull(partition_assignment) %>% 
  unique() %>% 
  sort()

# marker_plot_list<-lapply(
#   X = seq_along(titles),
#   FUN = plot_cells_alt,
#   cds = cds_aligned_all,
#   gene_or_genes = marker_gene_list,
#   h = 2.5,
#   w = 7.5,
#   cell_size = 0.5,
#   alpha = 0.5 ,
#   outfile = NULL,
#   plot_title = titles,
#   ncol = NULL,
#   plot_type = "png"
# )

gene_dot_plot1<-plot_genes_by_group(cds_aligned_all,
                                    unlist(marker_gene_list),
                                    group_cells_by="partition_assignment",
                                    ordering_type = "cluster_row_col",
                                    max.size=3,) + labs(x = "Partition Assignment")

save_plot(gene_dot_plot1, filename = "plots_out/gene_dot_plot1.pdf",base_width = 5, base_height = 5)

gene_dot_plot2<-plot_genes_by_group(cds_aligned_all,
                    top_specific_marker_ids,
                    group_cells_by="partition_assignment",
                    ordering_type="cluster_row_col",
                    max.size=3,) + labs(x = "Partition Assignment")
save_plot(gene_dot_plot2, filename = "plots_out/gene_dot_plot2.pdf", base_width = 5, base_height = 5)

three_gene<-plot_cells_alt(cds_aligned_all, gene_or_genes = c("CD34","MPO","FLT3"))
save_plot(three_gene, filename = "plots_out/three_gene.pdf", base_width = 7.5, base_height = 2.5)

# general group stats
## cell counts by cluster
granular_counts <- tbl_df(colData(cds_aligned_all)) %>%
  group_by(sample_id,partition_assignment) %>%
  summarise(
    patient = first(patient),
    timepoint = first(timepoint),
    tissue = first(tissue),
    n_cells = n()
  ) %>%
  ungroup() %>%
  mutate(running_total = cumsum(n_cells)) %>%
  write_csv(path = "data_out/partition_assignment_counts.csv")

## cell counts by sample
patient_timepoint_counts <- tbl_df(colData(cds_aligned_all)) %>%
  group_by(sample_id) %>%
  summarise(
    patient = first(patient),
    timepoint = first(timepoint),
    tissue = first(tissue),
    n_cells = n()
  ) %>%
  ungroup() %>%
  mutate(running_total = cumsum(n_cells)) %>%
  write_csv(path = "data_out/patient_timepoint_counts.csv")

save.image.pigz(file = "flt3_aml_bakerlab.RData", n.cores = 39)
