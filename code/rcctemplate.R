# Biobase is for accessing expressionSet objects
library(Biobase)

# Load Buettner data
devtools::install_github("jhsiao999/singleCellRNASeqHumanTungiPSC")
library(singleCellRNASeqHumanTungiPSC)

# Extract expression data
eset <- get(data("HumanTungiPSC"))

# Access the count matrix
counts <- exprs(eset)

# to get the logcount
# log_counts = log(counts)
# log_counts[which(log_counts == -Inf)] = NA
library(limma)
library(ashr)
source("~/HG/flash/Rcode/flash.R")
source("~/HG/flash/Rcode/greedy.R")
source("~/HG/flash/Rcode/backfitting.R")
voom_ouptut <- voom(counts, normalization = "none")
cpm_counts = voom_ouptut$E
log_cpm = log(cpm_counts)

g_flash = flashr::greedy(log_cpm, K = 30, flash_para = list(partype = "var_col"))
saveRDS(g_flash, "./logcpm/gflash_var_col.rds")
b_flash = flashr::backfitting(Y,initial_list = g_flash, flash_para = list(partype = "var_col"), maxiter_bf = 6)
saveRDS(b_flash, "./logcpm/bflash_var_col.rds")





source("~/HG/flash/Rcode/flash.R")
source("~/HG/flash/Rcode/greedy.R")
source("~/HG/flash/Rcode/backfitting.R")

saveRDS(result, "./NBSFAout/output.rds")


#!/bin/bash

#SBATCH --job-name=arrayJob
#SBATCH --output=/home/weidong/HG/flash/data/singlecell/HumanTungiPSC/outlog/arrayJob_%A_%a.out
#SBATCH --error=/home/weidong/HG/flash/data/singlecell/HumanTungiPSC/outlog/arrayJob_%A_%a.err
#SBATCH --array=1
#SBATCH --time=32:00:00
#SBATCH --partition=broadwl
#SBATCH --ntasks=1
#SBATCH --mem-per-cpu=4000


######################
# Begin work section #
######################

# Print this sub-job's task ID
cd /home/weidong/HG/flash/data/singlecell/HumanTungiPSC
Rscript --verbose /home/weidong/HG/flash/data/singlecell/HumanTungiPSC/logcpm.R

