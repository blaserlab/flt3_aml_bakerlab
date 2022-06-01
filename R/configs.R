# conflicts ---------------------------------------------------------------
# resolve conflicting function names here

conflict_prefer("filter", "dplyr")
conflict_prefer("lag", "dplyr")
conflict_prefer("rename", "dplyr")
conflict_prefer("select", "dplyr")
conflict_prefer("count", "dplyr")

theme_set(theme_cowplot(font_size = 10))

# analysis configurations -------------------------------------------------

experimental_group_palette = c(
  "Marrow T" = "#CCEBC5",
  "CD14 Monocytoid" = "#FFED6F",
  "Myeloblast" = "#BEBADA",
  "DC-like" = "#FB8072",
  "PC" = "#80B1D3",
  "B" = "#FDB462",
  "CD16 Monocytoid" = "#B3DE69",
  "Peripheral T" = "#FCCDE5",
  "PDC" = "#D9D9D9", 
  "WT" = "#000000", 
  "KO12" = "#0000FF", 
  "KO13" = "#A00000",
  "DMSO" = "#000000", 
  "CHMFL" = "#0000FF",
  "2" = brewer.pal(n = 7, name = "Set1")[1],
  "4" = brewer.pal(n = 7, name = "Set1")[2],
  "5" = brewer.pal(n = 7, name = "Set1")[3],
  "17" = brewer.pal(n = 7, name = "Set1")[4],
  "20" = brewer.pal(n = 7, name = "Set1")[5],
  "21" = brewer.pal(n = 7, name = "Set1")[6],
  "22" = brewer.pal(n = 7, name = "Set1")[7],
  "bm" = "#3C5488", 
  "blood" = "#DC0000",
  "sensitive" = "#7FBC41", 
  "unresponsive" = "#F1B6DA",
  "pre-treatment" = "grey90",
  "post-treatment" = "grey40"
)



network_out <- "/home/OSUMC.EDU/blas02/network/T/Labs/Baker-Sparreboom/FLT3-ITD Patient Samples/Gilteritinib/Manuscript/figs/source"
tables_directory <- "/home/OSUMC.EDU/blas02/network/T/Labs/Baker-Sparreboom/FLT3-ITD Patient Samples/Gilteritinib/Manuscript/tables"


# source local configs ----------------------------------------------------
# these are sourced after main configs and will overwrite duplicate entries if
# present. The file local_configs.R is ignored by git and so is useful for user-
# specific configurations such as output directories or formatting.

fs::file_create("R/local_configs.R") # will not overwrite

source("R/local_configs.R")
