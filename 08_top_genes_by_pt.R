source('00_packages_functions.R', echo=TRUE)

pt_list <- unique(colData(cds_aligned_all)$patient)
pt_list<-factor(pt_list, levels = c("1","2","3","4"))

colData(cds_aligned_all)$timepoint_binary<-recode(colData(cds_aligned_all)$timepoint, 
                                                  "C1D26_gilt" = "post-treatment",
                                                  "post_2_cycles_gilt" = "post-treatment",
                                                  "relapse" = "post-treatment",
                                                  "post_1_cycle_gilt" = "post-treatment",
                                                  "pretreatment" = "pre-treatment")

cds_no_lymphs<-cds_aligned_all[,partitions(cds_aligned_all) %notin% c(1,3,4,6)]

subset_get_markers <-
  function(cds,
           column_to_subset,
           value_to_select,
           group_var,
           outfile,
           n_genes,
           i = NULL) {
    if (!is.null(i)) {
      outfile<-outfile[i]
    } 
    cds_subset <- cds[,colData(cds)[[column_to_subset]] == value_to_select]
    tm <-
      top_markers(
        cds = cds_subset,
        group_cells_by = group_var,
        genes_to_test_per_group = n_genes,
        reference_cells = 1000,
        cores = 39
      )
    tm %>% arrange(cell_group) %>% write_csv(path = outfile)
  }

class(colData(cds_no_lymphs)[["patient"]])
pt_list

lapply(X = seq_along(pt_list), 
       FUN = subset_get_markers,
       cds = cds_no_lymphs,
       column_to_subset = "patient",
       value_to_select = pt_list,
       group_var = "timepoint_binary",
       outfile = paste0("data_out/topgenes_prepost/pt_",pt_list,".csv"),
       n_genes = 50)


save.image.pigz(file = "flt3_aml_bakerlab.RData", n.cores = 39)
