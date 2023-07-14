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


# if using the blaseRtemplates environment, edit this line and run to load the data package-------------------------------------
blaseRtemplates::project_data(path = "/home/OSUMC.EDU/blas02/network/X/Labs/Blaser/share/collaborators/bakerlab_flt3_gilt_project/datapkg")

# if using standard R data functions, uncomment and run these line to load the data package into your global environment:


d <- data(package = "flt3.aml.bakerlab.datapkg")
data(list = d$results[, "Item"], package = "flt3.aml.bakerlab.datapkg")
rm(d)



