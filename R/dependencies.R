# renv --------------------------------------------------------------------

# set up the renv from scratch

# renv::init(bioconductor = TRUE)

# restore the renv from the lockfile

# renv::restore()



# package installation ----------------------------------------------------

# # Try this first...it's faster:
# blaseRtemplates::easy_install("CellChat", how = "link_from_cache")

# # If you need a new package or an update, try this:
# blaseRtemplates::easy_install("", how = "new_or_update")
# blaseRtemplates::easy_install("blaserlab/blaseRtools", how = "new_or_update")

# # If you are installing from a "tarball", use this:
# blaseRtemplates::easy_install("/path/to/tarball.tar.gz")

# # use "bioc::<package name>" for bioconductor packages
# # use "<repo/package name>" for github source packages

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

# uncomment and use the following to install or update the data package---------------------------------------
bb_renv_datapkg(path = "/home/OSUMC.EDU/blas02/network/X/Labs/Blaser/share/collaborators/bakerlab_flt3_gilt_project/datapkg")

# use this to load the data package-------------------------------------

# run once to load, run again to unload
requireData(package = "flt3.aml.bakerlab.datapkg")
