source("R/configs.R")

pseudobulk_pre_myeloblast %>%
  mutate(up_in = ifelse(log2FoldChange>0,"unresponsive","sensitive")) %>%
  relocate(gene_short_name, up_in) %>% #filter(padj<0.05) %>% View
  arrange(padj) %>%
  write_tsv(str_glue("{tables_directory}/myeloblasts_pre_treatment.tsv"))

pseudobulk_post_myeloblast %>%
  mutate(up_in = ifelse(log2FoldChange>0,"unresponsive","sensitive")) %>%
  relocate(gene_short_name, up_in) %>% #filter(padj<0.05) %>% View
  arrange(padj) %>% 
  write_tsv(str_glue("{tables_directory}/myeloblasts_post_treatment.tsv"))

partition_tm_res_cds_aligned %>%
  as_tibble() %>%
  left_join(
    colData(cds_anno_aligned_tissue_id) %>%
      as_tibble() %>%
      group_by(partition_assignment, partition_assignment_1) %>%
      summarise()
    ) %>%
  relocate(partition_assignment_1) %>%
  select(-c(partition_assignment, cell_group)) %>%
  arrange(partition_assignment_1, marker_test_q_value) %>%
  dplyr::rename(partition = partition_assignment_1) %>%
  write_tsv(str_glue("{tables_directory}/partition_top_markers.tsv"))

rowData(cds_anno_aligned_tissue_id) %>%
  as_tibble() %>%
  select(id, gene_short_name, module) %>%
  arrange(module, gene_short_name) %>%
  # write_tsv(str_glue("{tables_directory}/gene_modules.tsv"))
  write_tsv("~/network/P/blaser_lab_p/writing/flt3_aml_bakerlab_manuscript/tables/table_s4.tsv")

gorilla_modules_cds_anno %>%
  arrange(module, `FDR q-value`) %>%
  select(module, `GO Term`, Description, `P-value`, `FDR q-value`, Enrichment) %>%
  write_tsv(str_glue("{tables_directory}/module_goterms.tsv"))
