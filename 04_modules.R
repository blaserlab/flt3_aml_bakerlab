source('00_packages_functions.R', echo=TRUE)
load.pigz(file = "flt3_aml_bakerlab.RData")

pr_graph_test_res <- graph_test(cds_aligned_all, neighbor_graph="knn", cores=39)
pr_deg_ids <- row.names(subset(pr_graph_test_res, q_value < 0.05))

gene_module_df <-
  find_gene_modules(cds_aligned_all[pr_deg_ids, ], resolution = 1e-2)

gene_module_anno<-left_join(gene_module_df,tbl_df(rowData(cds_aligned_all)))

View(gene_module_anno %>% filter(module==20))

colData(cds_aligned_all)$pa_r <-
  paste0(
    colData(cds_aligned_all)$partition_assignment,
    " ",
    colData(cds_aligned_all)$response
  )

colData(cds_aligned_all)$paptt <-
  paste0(
    colData(cds_aligned_all)$partition_assignment,
    " ",
    colData(cds_aligned_all)$pt_response,
    " ",
    colData(cds_aligned_all)$timepoint
  )
cell_group_df <-
  tibble::tibble(cell = row.names(colData(cds_aligned_all)),
                 cell_group = colData(cds_aligned_all)$paptt)

agg_mat <- aggregate_gene_expression(cds_aligned_all, gene_module_df, cell_group_df)
row.names(agg_mat) <- stringr::str_c("Module ", row.names(agg_mat))
colnames(agg_mat) <- stringr::str_c(colnames(agg_mat))

modules_granular<-pheatmap::pheatmap(agg_mat, cluster_rows=TRUE, cluster_cols=TRUE,
                   scale="column", clustering_method="ward.D2",
                   fontsize=10)
# redo modules excluding B1, PC, T, T/NK
agg_mat1 <- aggregate_gene_expression(cds_aligned_all[,colData(cds_aligned_all)$partition %notin% c(1,3,4,6)], gene_module_df, cell_group_df)
row.names(agg_mat1) <- stringr::str_c("Module ", row.names(agg_mat1))
colnames(agg_mat1) <- stringr::str_c(colnames(agg_mat1))

#set up the annotation dataframe
annotation_col0 = tibble( # this makes a dataframe with a line for each row. The annotation variable will be "OfNote"
  columns = colnames(agg_mat1)# this puts the columns from the module plot into rows
)

annotation_col <- left_join(
  annotation_col0,
  tbl_df(colData(cds_aligned_all)) %>%
    select(partition_assignment, paptt),
  by = c("columns" = "paptt")
) %>%
  unique()%>%
  rename("Partition" = "partition_assignment") %>%
  column_to_rownames(var = "columns")

Partition<-brewer.pal(n = 7, name = "Set1")
names(Partition) <- annotation_col %>% pull(Partition) %>% unique()
anno_colors <- list(Partition = Partition)
modules_granular_noTBNKPC <- as.ggplot(
  pheatmap::pheatmap(
    agg_mat1,
    cluster_rows = TRUE,
    cluster_cols = TRUE,
    scale = "column",
    clustering_method = "ward.D2",
    fontsize = 6,
    annotation_col = annotation_col,
    annotation_colors = anno_colors
    
  )
)
save_plot(modules_granular_noTBNKPC, file = "plots_out/modules_granular_noTBNKPC.pdf", base_height = 10, base_width = 7.5)

#redo the module heatmap with jut aml/prog and blast1

cds_aligned_aml_prog_blast1<-cds_aligned_all[,colData(cds_aligned_all)$partition_assignment %in% c("AML/Progenitor","Blast 1")]
cell_group_df_aml_prog_blast1 <-
  tibble::tibble(cell = row.names(colData(cds_aligned_aml_prog_blast1)),
                 cell_group = colData(cds_aligned_aml_prog_blast1)$paptt)

agg_mat_aml_prog_blast1 <- aggregate_gene_expression(cds_aligned_aml_prog_blast1, gene_module_df, cell_group_df_aml_prog_blast1)
row.names(agg_mat_aml_prog_blast1) <- stringr::str_c("Module ", row.names(agg_mat_aml_prog_blast1))
colnames(agg_mat_aml_prog_blast1) <- stringr::str_c(colnames(agg_mat_aml_prog_blast1))

annotation_col0 = tibble( # this makes a dataframe with a line for each row. The annotation variable will be "OfNote"
  columns = colnames(agg_mat_aml_prog_blast1)# this puts the columns from the module plot into rows
)

annotation_col_aml_prog_blast1 <- left_join(
  annotation_col0,
  tbl_df(colData(cds_aligned_aml_prog_blast1)) %>%
    select(partition_assignment, paptt,response),
  by = c("columns" = "paptt")
) %>%
  unique()%>%
  rename("Partition" = "partition_assignment","Response" = "response") %>%
  column_to_rownames(var = "columns")

Partition<-c("#E41A1C","#377EB8")
names(Partition) <- annotation_col_aml_prog_blast1 %>% pull(Partition) %>% unique()
Response<-c("#4DAF4A","#984EA3","#FF7F00")
names(Response)<-annotation_col_aml_prog_blast1 %>% pull(Response) %>% unique() %>% sort()
anno_colors <- list(Partition = Partition,Response = Response)

modules_granular_noTBNKPC <- as.ggplot(
  pheatmap::pheatmap(
    agg_mat_aml_prog_blast1,
    cluster_rows = TRUE,
    cluster_cols = TRUE,
    scale = "column",
    clustering_method = "ward.D2",
    fontsize = 6,
    annotation_col = annotation_col_aml_prog_blast1,
    annotation_colors = anno_colors
    
  )
)
save_plot(modules_granular_noTBNKPC, file = "plots_out/modules_aml_blast1.pdf", base_height = 10, base_width = 7.5)


gene_module_anno %>% write_csv(path = "data_out/gene_modules.csv")


save.image.pigz(file = "flt3_aml_bakerlab.RData", n.cores = 39)
