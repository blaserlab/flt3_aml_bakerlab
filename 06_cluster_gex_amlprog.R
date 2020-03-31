source('00_packages_functions.R', echo=TRUE)
#load.pigz(file = "flt3_aml_bakerlab.RData")

# select amlprog clusters and plot
cds_amlprog<-cds_aligned_all[,colData(cds_aligned_all)$partition_assignment=="AML/Progenitor"]

#focusing on original clusters
amlprog_clusters <-
  custom_cp_plot(
    cds_amlprog,
    cp = "cluster",
    group_label_size = 4,
    plot_title = "",
    outfile = "plots_out/amlprog_clusters.pdf",
    w = 3.3,
    h = 3
  )

# stratify by patient
titles<-unique(colData(cds_amlprog)$pt_response)
titles<-factor(titles, levels = c("Pt 1:  PR","Pt 2:  refractory", "Pt 3:  CR", "Pt 4:  refractory"))
titles<-sort(titles)
outfiles<-c("amlprog_pt1_pr","amlprog_pt2_ref","amlprog_pt3_cr","amlprog_pt4_ref")

amlprog_by_pt_list <-
  lapply(
    X = seq_along(titles),
    FUN = custom_variable_plot,
    cds = cds_amlprog,
    var = "pt_response",
    palette_viridis = F,
    value_to_highlight = titles,
    legend_pos = "none",
    cell_size = 1,
    outfile = outfiles,
    h = 3,
    w = 3.3,
    plot_title = titles,
    foreground_alpha = 0.2,
    legend_title = NULL
  )


# generate new categorical variable for cluster 71 and not 71

colData(cds_amlprog)$cluster_trinary <-
  recode(
    colData(cds_amlprog)$cluster,
    "30" = "pt1_pt4",
    "43" = "pt1_pt4",
    "62" = "pt1_pt4",
    "59" = "pt1_pt4",
    "58" = "pt1_pt4",
    "57" = "pt1_pt4",
    "56" = "pt1_pt4",
    "55" = "pt1_pt4",
    "54" = "pt1_pt4",
    "53" = "pt1_pt4",
    "70" = "pt1_pt4",
    "6" = "pt3",
    "31" = "pt3",
    "32" = "pt3",
    "33" = "pt3",
    "34" = "pt3",
    "35" = "pt3",
    "42" = "pt3",
    "47" = "pt3",
    "5" = "all",
    "29" = "all",
    "44" = "all",
    "27" = "all",
    "74" = "all",
    "67" = "all",
    "8" = "all",
    "4" = "all",
    "20" = "all",
    "2" = "all",
    "28" = "all"
  )

colData(cds_amlprog)$cluster_trinary <-
  factor(colData(cds_amlprog)$cluster_trinary,
         levels = c("pt1_pt4", "pt3", "all"))
amlprog_clusters_trinary <-
  custom_variable_plot(
    cds_amlprog,
    var = "cluster_trinary",
    outfile = "plots_out/amlprog_clusters_trinary.pdf",
    foreground_alpha = 0.2,
    w = 3.8,
    h = 3
  )


#now identify top markers based on clusters trinary
amlprog_trinary_tm<-top_markers(
  cds = cds_amlprog,
  group_cells_by = "cluster_trinary",
  genes_to_test_per_group = 25,
  reference_cells = 1000,
  cores = 39,
  verbose = T
)

pt1_pt4_genes<-amlprog_trinary_tm %>% filter(cell_group == "pt1_pt4") %>% pull(gene_short_name)

pt3_genes<-amlprog_trinary_tm %>% filter(cell_group == "pt3") %>% pull(gene_short_name)

all_genes<-amlprog_trinary_tm %>% filter(cell_group == "all") %>% pull(gene_short_name)


#make violins
violin_up_in_pt1_pt4 <- custom_violin_plot(
  cds = cds_amlprog, 
  plot_title = "Genes up in Pt 1 (PR) and Pt 4 (Refractory)",
  genes_to_plot = pt1_pt4_genes,
  variable = "cluster_trinary",
  rows = 5,
  outfile = "plots_out/violin_up_in_pt1_pt4.pdf",
  w = 7.5,
  h = 8.5,
)

violin_up_in_pt3 <- custom_violin_plot(
  cds = cds_amlprog, 
  plot_title = "Genes up in Pt 3 (CR)",
  genes_to_plot = pt3_genes,
  variable = "cluster_trinary",
  rows = 5,
  outfile = "plots_out/violin_up_in_pt3.pdf",
  w = 7.5,
  h = 8.5,
)

violin_up_in_all_amlprog <- custom_violin_plot(
  cds = cds_amlprog, 
  plot_title = "Genes up in All",
  genes_to_plot = all_genes,
  variable = "cluster_trinary",
  rows = 5,
  outfile = "plots_out/violin_up_in_all.pdf",
  w = 7.5,
  h = 8.5,
)


#save.image.pigz(file = "flt3_aml_bakerlab.RData", n.cores = 39)
