source('00_packages_functions.R', echo=TRUE)

#load in patient metadata
md0 <- read_csv("scRNASeq_samples.csv") 
md <- md0 %>%
  group_by(patient) %>%
  summarise(day_0 = min(sample_date)) %>%
  left_join(md0,.) %>%
  mutate(treatment_day = (sample_date - day_0))
md

# load in 10X pipestance:  do this by :pwd in vifm, then highlight the file path.  Will be copied to clipboard.
cds_Y5208 <-
  load_cellranger_data(pipestance_path = "/home/OSUMC.EDU/blas02/network/X/Labs/Blaser/single_cell/baker_lab_dec2019/output_baker_dec2019/baker_dec2019_Y5208_BakerS_2305_V1G_1", barcode_filtered = TRUE)
colData(cds_Y5208)$sample_id<-"Y5208"
cds_Y5209 <-
  load_cellranger_data(pipestance_path = "/home/OSUMC.EDU/blas02/network/X/Labs/Blaser/single_cell/baker_lab_dec2019/output_baker_dec2019/baker_dec2019_Y5209_BakerS_4592_V1G_1", barcode_filtered = TRUE)
colData(cds_Y5209)$sample_id<-"Y5209"
cds_Y6584 <-
  load_cellranger_data(pipestance_path = "/home/OSUMC.EDU/blas02/network/X/Labs/Blaser/single_cell/baker_lab_mar2020/output_baker_1560_1561_mar2020/baker_1560_mar2020_Y6584_BakerS_1_V1G_1_GEX", barcode_filtered = TRUE)
colData(cds_Y6584)$sample_id<-"Y6584"
cds_Y6585 <-
  load_cellranger_data(pipestance_path = "/home/OSUMC.EDU/blas02/network/X/Labs/Blaser/single_cell/baker_lab_mar2020/output_baker_1560_1561_mar2020/baker_1560_mar2020_Y6585_BakerS_2_V1G_1_GEX", barcode_filtered = TRUE)
colData(cds_Y6585)$sample_id<-"Y6585"
cds_Y5206 <-
  load_cellranger_data(pipestance_path = "/home/OSUMC.EDU/blas02/network/X/Labs/Blaser/single_cell/baker_lab_dec2019/output_baker_dec2019/baker_dec2019_Y5206_BakerS_1110_V1G_1", barcode_filtered = TRUE)
colData(cds_Y5206)$sample_id<-"Y5206"
cds_Y5207 <-
  load_cellranger_data(pipestance_path = "/home/OSUMC.EDU/blas02/network/X/Labs/Blaser/single_cell/baker_lab_dec2019/output_baker_dec2019/baker_dec2019_Y5207_BakerS_1669_V1G_1", barcode_filtered = TRUE)
colData(cds_Y5207)$sample_id<-"Y5207"
cds_Y6586 <-
  load_cellranger_data(pipestance_path = "/home/OSUMC.EDU/blas02/network/X/Labs/Blaser/single_cell/baker_lab_mar2020/output_baker_1560_1561_mar2020/baker_1561_mar2020_Y6586_BakerS_3_V1G_1_GEX", barcode_filtered = TRUE)
colData(cds_Y6586)$sample_id<-"Y6586"
cds_Y6587 <-
  load_cellranger_data(pipestance_path = "/home/OSUMC.EDU/blas02/network/X/Labs/Blaser/single_cell/baker_lab_mar2020/output_baker_1560_1561_mar2020/baker_1561_mar2020_Y6587_BakerS_4_V1G_1_GEX", barcode_filtered = TRUE)
colData(cds_Y6587)$sample_id<-"Y6587"
cds_Y6588 <-
  load_cellranger_data(pipestance_path = "/home/OSUMC.EDU/blas02/network/X/Labs/Blaser/single_cell/baker_lab_mar2020/output_baker_1560_1561_mar2020/baker_1561_mar2020_Y6588_BakerS_5_V1G_1_GEX", barcode_filtered = TRUE)
colData(cds_Y6588)$sample_id<-"Y6588"
cds_Y6589 <-
  load_cellranger_data(pipestance_path = "/home/OSUMC.EDU/blas02/network/X/Labs/Blaser/single_cell/baker_lab_mar2020/output_baker_1560_1561_mar2020/baker_1561_mar2020_Y6589_BakerS_6_V1G_1_GEX", barcode_filtered = TRUE)
colData(cds_Y6589)$sample_id<-"Y6589"


# generate list of cds to combine and then combine
cds_list_all <-
  list(
    cds_Y5206,
    cds_Y5207,
    cds_Y5208,
    cds_Y5209,
    cds_Y6584,
    cds_Y6585,
    cds_Y6586,
    cds_Y6587,
    cds_Y6588,
    cds_Y6589
  )
cds_list_onlybm <-
  list(
    cds_Y5206,
    cds_Y5207,
    cds_Y5208,
    cds_Y5209,
    cds_Y6586,
    cds_Y6587,
    cds_Y6588,
    cds_Y6589
  )

cds_combined_all <-
  combine_cds(cds_list = cds_list_all, keep_all_genes = TRUE)
cds_combined_onlybm <-
  combine_cds(cds_list = cds_list_onlybm, keep_all_genes = TRUE)

# join the metadata onto the cds's
colData(cds_combined_all)[5:14]<-left_join(tbl_df(colData(cds_combined_all)),md) %>% select(c(5:14)) 
colData(cds_combined_onlybm)[5:14]<-left_join(tbl_df(colData(cds_combined_onlybm)),md) %>% select(c(5:14)) 


#optional - trim off uninformative genes
cds_trimmed_all <-
  cds_combined_all[substr(rowData(cds_combined_all)$gene_short_name, 1, 2) !=
                     "RP", ]
cds_trimmed_onlybm <-
  cds_combined_onlybm[substr(rowData(cds_combined_onlybm)$gene_short_name, 1, 2) !=
                        "RP", ]

## Pre-process the data
cds_trimmed_all <- preprocess_cds(cds_trimmed_all, num_dim = 100,)
cds_trimmed_onlybm <-
  preprocess_cds(cds_trimmed_onlybm, num_dim = 100)

# batch correct based on the patient variable
cds_aligned_all<-align_cds(cds_trimmed_all, alignment_group = "patient")
cds_aligned_onlybm<-align_cds(cds_trimmed_onlybm, alignment_group = "patient")

## Reduce dimensionality
cds_trimmed_all<-reduce_dimension(cds_trimmed_all, cores = 39)
cds_trimmed_onlybm<-reduce_dimension(cds_trimmed_onlybm, cores = 39)

cds_aligned_all<-reduce_dimension(cds_aligned_all, cores = 39)
cds_aligned_onlybm<-reduce_dimension(cds_aligned_onlybm, cores = 39)

save.image.pigz(file = "flt3_aml_bakerlab.RData", n.cores = 39)
