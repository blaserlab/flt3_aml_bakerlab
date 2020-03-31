source('00_packages_functions.R', echo=TRUE)
#load.pigz(file = "flt3_aml_bakerlab.RData")

blast1_binary_top100<-top_markers(
  cds = cds_blast1,
  group_cells_by = "cluster_binary",
  genes_to_test_per_group = 100,
  reference_cells = 1000,
  cores = 39,
  verbose = T
)

blast1_binary_top100 %>% tbl_df() %>% arrange(cell_group) %>% write_csv("data_out/blast1_binary_top100.csv")

amlprog_trinary_top100<-top_markers(
  cds = cds_amlprog,
  group_cells_by = "cluster_trinary",
  genes_to_test_per_group = 100,
  reference_cells = 1000,
  cores = 39,
  verbose = T
)

amlprog_trinary_top100 %>% tbl_df() %>% arrange(cell_group) %>% write_csv("data_out/amlprog_trinary_top100.csv")

save.image.pigz(file = "flt3_aml_bakerlab.RData", n.cores = 39)
