source('00_packages_functions.R', echo=TRUE)
#load.pigz(file = "flt3_aml_bakerlab.RData")

# select blast 1 clusters and plot
cds_blast1<-cds_aligned_all[,colData(cds_aligned_all)$partition_assignment=="Blast 1"]

#focusing on original clusters
blast1_clusters <-
  custom_cp_plot(
    cds_blast1,
    cp = "cluster",
    group_label_size = 4,
    plot_title = "",
    outfile = "plots_out/blast1_clusters.pdf",
    w = 3.3,
    h = 3
  )

# focus on patient 3
pt3_blast1 <- custom_variable_plot(
  cds_blast1,
  var = "patient",
  palette_viridis = F,
  value_to_highlight = 3,
  legend_pos = "none",
  cell_size = 1,
  outfile = "plots_out/pt3_blast1.pdf",
  h = 3,
  w = 3.3,
  plot_title = "Patient 3"
)

# generate new categorical variable for cluster 71 and not 71

colData(cds_blast1)$cluster_binary <-
  recode(
    colData(cds_blast1)$cluster,
    "71" = "71",
    "23" = "not71",
    "24" = "not71",
    "25" = "not71",
    "26" = "not71",
    "61" = "not71",
    "63" = "not71",
    "64" = "not71",
    "65" = "not71",
    "68" = "not71",
    "72" = "not71",
    "73" = "not71",
    "75" = "not71"
  )

blast1_clusters_binary <-
  custom_variable_plot(
    cds_blast1,
    var = "cluster_binary",
    outfile = "plots_out/blast1_clusters_binary.pdf",
    w = 3.8,
    h = 3
  )


#now identify top markers based on clusters binary
blast1_binary_tm<-top_markers(
  cds = cds_blast1,
  group_cells_by = "cluster_binary",
  genes_to_test_per_group = 25,
  reference_cells = 1000,
  cores = 39,
  verbose = T
)

not71genes<-blast1_binary_tm %>% filter(cell_group == "not71") %>% pull(gene_short_name)

yes71genes<-blast1_binary_tm %>% filter(cell_group == "71") %>% pull(gene_short_name)

#make violins
violin_up_in_71 <- custom_violin_plot(
  cds = cds_blast1,
  genes_to_plot = yes71genes,
  variable = "cluster_binary",
  plot_title = "Genes up in Cluster 71",
  rows = 5,
  outfile = "plots_out/violin_up_in_71.pdf",
  w = 5,
  h = 8.5
)

violin_up_in_not71 <- custom_violin_plot(
  cds = cds_blast1,
  genes_to_plot = not71genes,
  variable = "cluster_binary",
  rows = 5,
  plot_title = "Genes up in Not Cluster 71",
  outfile = "plots_out/violin_up_in_not71.pdf",
  w = 5,
  h = 8.5
)

#save.image.pigz(file = "flt3_aml_bakerlab.RData", n.cores = 39)
