#!/usr/bin/Rscript

# --------------------------------------------------------------------------
# Fetch human fetal data and generate human chondrocyte signature
# --------------------------------------------------------------------------

# Set adaptive working directory
if(dir.exists("/Users/mac/cloudstor/")) {
  setwd("/Users/mac/cloudstor/scott_projects/osteochondro_signature/") # Uses practice data (5% of cells from each sample) if running locally
  place <- "local"
} else {
  setwd("/share/ScratchGeneral/scoyou/sarah_projects/osteochondro_signature/")
  place <- "wolfpack"
}

library('ExpressionAtlas')


# Load counts data for human tissues from array express (note: GTEx was used to supplement fetal tissues)
# --------------------------------------------------------------------------

# GTEx data of human tissues
lookup <- data.frame(searchAtlasExperiments( c('GTEx'), species = " Homo sapiens" ))
gtex <- getAtlasData("E-MTAB-5214")
sumexp <- gtex[[ "E-MTAB-5214" ]]$rnaseq
gtex_counts <- assays( sumexp )$counts

# Fetal tissues
lookup <- data.frame(searchAtlasExperiments( c('fetus'), species = " Homo sapiens" ))
data <- getAtlasData(c("E-MTAB-3871", "E-GEOD-59089", "E-MTAB-4840"))
mtab3871 <- data[[ "E-MTAB-3871" ]]$rnaseq
geod59089 <- data[[ "E-GEOD-59089" ]]$rnaseq
mtab4840 <- data[[ "E-MTAB-4840" ]]$rnaseq





