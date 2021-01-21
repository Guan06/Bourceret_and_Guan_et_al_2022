#!/biodata/dep_psl/grp_rgo/tools/bin/Rscript

library(ggplot2)
library(parallelDist)
source("../00.common_scripts/plot_settings.R")
source("../00.common_scripts/pcoa_plotting.R")

###############################################################################
###############################################################################
## For Bacteria

asv_file <- "../00.data/final_ASV/Bacteria/ASV_raref.rds"
design_file <- "../00.data/design.txt"
kingdom <- "Bacteria"

asv <- readRDS(asv_file)
design <- read.table(design_file, header = T, sep = "\t")

# Get the intersection of asv and design file
inter <- intersect(as.character(colnames(asv)), as.character(design$Sample_ID))
asv <- asv[, colnames(asv) %in% inter]
design <- design[design$Sample_ID %in% inter, ]

# Bray-Curtis distance
bc <- as.matrix(parDist(t(asv), method = "bray"))
des <- design
dis <- bc
dmr <- cmdscale(dis, k = 4, eig = T)

p1 <- pcoa(dmr, des, 12, "Management", "Stage", 1.2, kingdom) +
    theme(legend.position = "none")

##############################################################################
###############################################################################
## For Fungi

asv_file <- "../00.data/final_ASV/Fungi/ASV_raref.rds"
kingdom <- "Fungi"

asv <- readRDS(asv_file)
design <- read.table(design_file, header = T, sep = "\t")

# Get the intersection of asv and design file
inter <- intersect(as.character(colnames(asv)), as.character(design$Sample_ID))
asv <- asv[, colnames(asv) %in% inter]
design <- design[design$Sample_ID %in% inter, ]

# Bray-Curtis distance
bc <- as.matrix(parDist(t(asv), method = "bray"))
des <- design
dis <- bc
dmr <- cmdscale(dis, k = 4, eig = T)

fp1 <- pcoa(dmr, des, 12, "Management", "Stage", 1.2, kingdom) +
    theme(legend.position = "none")

###############################################################################
###############################################################################
## For oomycetes

asv_file <- "../00.data/final_ASV/Oomycetes/ASV_raref.rds"
kingdom <- "Oomycetes"

asv <- readRDS(asv_file)
design <- read.table(design_file, header = T, sep = "\t")

# Get the intersection of asv and design file
inter <- intersect(as.character(colnames(asv)), as.character(design$Sample_ID))
asv <- asv[, colnames(asv) %in% inter]
design <- design[design$Sample_ID %in% inter, ]

# Bray-Curtis distance
bc <- as.matrix(parDist(t(asv), method = "bray"))
des <- design
dis <- bc
dmr <- cmdscale(dis, k = 4, eig = T)
op1 <- pcoa(dmr, des, 12, "Management", "Stage", 1.2, kingdom) +
    theme(legend.position = "none")

###############################################################################
## Output the files
library(gridExtra)

all_p123 <- grid.arrange(p1, fp1, op1, nrow = 1, ncol = 3)
ggsave("Figure_S3.pdf", all_p123, width = 10, height = 3.33)
