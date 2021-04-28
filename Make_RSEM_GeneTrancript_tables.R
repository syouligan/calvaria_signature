#!/usr/bin/Rscript

# --------------------------------------------------------------------------
# Make tables of RSEM results
# --------------------------------------------------------------------------

# Set adaptive working directory
if(dir.exists("/Users/mac/cloudstor/")) {
  setwd("/Users/mac/cloudstor/scott_projects/osteochondro_signature/") # Uses practice data (5% of cells from each sample) if running locally
  place <- "local"
} else {
  setwd("/share/ScratchGeneral/scoyou/sarah_projects/osteochondro_signature/")
  place <- "wolfpack"
}

library('data.table')

#Setup in/out directory structure 
# -------------------------------

# In path
samplePath <- "/share/ScratchGeneral/scoyou/scott_projects/osteochondro_signature/data/humanchondro_bulk/project_results/RSEM"

# Outpath
outPath <- "/share/ScratchGeneral/scoyou/scott_projects/osteochondro_signature/data/humanchondro_bulk/project_results/RSEM_tables/"
dir.create("outPath")

# Make tables containing RSEM gene counts, TMP, FPKM per sample
# -------------------------------

# List files and sample names
sampleFiles <- list.files(path=samplePath, pattern= "genes\\.results$", recursive=TRUE, full.names=TRUE)
sampleNames <- list.files(path=samplePath)

# Make FPKM table
RSEM_FPKM <- data.frame(fread(sampleFiles[1], header=TRUE, select=c(1)), row.names=1)
for (i in sampleFiles) {
  RSEM_FPKM <- cbind(RSEM_FPKM, fread(i, header=TRUE, select=c(7)))
}
colnames(RSEM_FPKM) <- sampleNames
write.csv(RSEM_FPKM, file=paste(outPath, "RSEM_FPKM_genes.csv", sep=""))

# Make TPM table
RSEM_TPM <- data.frame(fread(sampleFiles[1], header=TRUE, select=c(1)), row.names=1)
for (i in sampleFiles) {
  RSEM_TPM <- cbind(RSEM_TPM, fread(i, header=TRUE, select=c(6)))
}
colnames(RSEM_TPM) <- sampleNames
write.csv(RSEM_TPM, file=paste(outPath, "RSEM_TPM_genes.csv", sep=""))

# Make counts table
RSEM_count <- data.frame(fread(sampleFiles[1], header=TRUE, select=c(1)), row.names=1)
for (i in sampleFiles) {
  RSEM_count <- cbind(RSEM_count, fread(i, header=TRUE, select=c(5)))
}
colnames(RSEM_count) <- sampleNames
RSEM_all_count <- RSEM_count
write.csv(RSEM_all_count, file=paste(outPath, "RSEM_all_count_genes.csv", sep=""))

# Make tables containing RSEM gene transcript, TMP, FPKM per sample
# -------------------------------

# List files and sample names
sampleFiles <- list.files(path=samplePath, pattern= "isoforms\\.results$", recursive=TRUE, full.names=TRUE)
sampleNames <- list.files(path=samplePath)

# Make FPKM table
RSEM_FPKM <- data.frame(fread(sampleFiles[1], header=TRUE, select=c(1)), row.names=1)
for (i in sampleFiles) {
  RSEM_FPKM <- cbind(RSEM_FPKM, fread(i, header=TRUE, select=c(7)))
}
colnames(RSEM_FPKM) <- sampleNames
write.csv(RSEM_FPKM, file=paste(outPath, "RSEM_FPKM_transcripts.csv", sep=""))

# Make TPM table
RSEM_TPM <- data.frame(fread(sampleFiles[1], header=TRUE, select=c(1)), row.names=1)
for (i in sampleFiles) {
  RSEM_TPM <- cbind(RSEM_TPM, fread(i, header=TRUE, select=c(6)))
}
colnames(RSEM_TPM) <- sampleNames
write.csv(RSEM_TPM, file=paste(outPath, "RSEM_TPM_transcripts.csv", sep=""))

# Make counts table
RSEM_count <- data.frame(fread(sampleFiles[1], header=TRUE, select=c(1)), row.names=1)
for (i in sampleFiles) {
  RSEM_count <- cbind(RSEM_count, fread(i, header=TRUE, select=c(5)))
}
colnames(RSEM_count) <- sampleNames
RSEM_all_count <- RSEM_count
write.csv(RSEM_all_count, file=paste(outPath, "RSEM_all_count_transcripts.csv", sep=""))

# Make IsoformPercent table
RSEM_IsoPct <- data.frame(fread(sampleFiles[1], header=TRUE, select=c(1)), row.names=1)
for (i in sampleFiles) {
  RSEM_IsoPct <- cbind(RSEM_IsoPct, fread(i, header=TRUE, select=c(8)))
}
colnames(RSEM_IsoPct) <- sampleNames
write.csv(RSEM_IsoPct, file=paste(outPath, "RSEM_IsoPct_transcripts.csv", sep=""))

