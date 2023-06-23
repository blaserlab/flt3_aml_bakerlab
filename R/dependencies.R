# load and attach packages ------------------------------------------------
suppressPackageStartupMessages(library("conflicted"))
suppressPackageStartupMessages(library("blaseRtools"))
suppressPackageStartupMessages(library("blaseRtemplates"))
suppressPackageStartupMessages(library("tidyverse"))
suppressPackageStartupMessages(library("monocle3"))
suppressPackageStartupMessages(library("circlize"))
suppressPackageStartupMessages(library("ComplexHeatmap"))
suppressPackageStartupMessages(library("lazyData"))
suppressPackageStartupMessages(library("cowplot"))
suppressPackageStartupMessages(library("RColorBrewer"))
suppressPackageStartupMessages(library("ggrepel"))
suppressPackageStartupMessages(library("ggpubr"))
suppressPackageStartupMessages(library("rstatix"))
suppressPackageStartupMessages(library("stats"))
suppressPackageStartupMessages(library("ggtext"))
suppressPackageStartupMessages(library("CellChat"))


# use this to load the data package-------------------------------------
blaseRtemplates::project_data(path = "/home/OSUMC.EDU/blas02/network/X/Labs/Blaser/share/collaborators/bakerlab_flt3_gilt_project/datapkg")

