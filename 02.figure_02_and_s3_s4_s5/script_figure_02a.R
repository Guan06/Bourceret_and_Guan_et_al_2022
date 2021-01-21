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

###############################################################################
# Plotting for A

des <- design
dis <- bc
dmr <- cmdscale(dis, k = 4, eig = T)
n <- nrow(dis)
print(n)

## Compartment + Stage / Management + Stage
p1 <- pcoa(dmr, des, 12, "Compartment", "Stage", 1.2, kingdom) +
    theme(legend.position = "none")

p2 <- pcoa(dmr, des, 13, "Compartment", "Stage", 1.2, kingdom) +
    theme(legend.position ="none")

p3 <- pcoa(dmr, des, 34, "Management", "Stage", 1.2, kingdom) +
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

###############################################################################
# Plotting for all samples
des <- design
dis <- bc
dmr <- cmdscale(dis, k = 4, eig = T)
n <- nrow(dis)
print(n)

## Comfpartment + Stage / Management + Stage
fp1 <- pcoa(dmr, des, 12, "Compartment", "Stage", 1.2, kingdom) +
    theme(legend.position = "none")

fp2 <- pcoa(dmr, des, 13, "Compartment", "Stage", 1.2, kingdom) +
    theme(legend.position = "none")

fp3 <- pcoa(dmr, des, 34, "Management", "Stage", 1.2, kingdom) +
    theme(legend.position = "none")

###############################################################################
###############################################################################
## For Oomycetes

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

###############################################################################
# Plotting for all samples
des <- design
dis <- bc
dmr <- cmdscale(dis, k = 4, eig = T)
n <- nrow(dis)
print(n)

## Compartment + Stage / Management + Stage
op1 <- pcoa(dmr, des, 12, "Compartment", "Stage", 1.2, kingdom) +
    theme(legend.position = "none")

op2 <- pcoa(dmr, des, 13, "Compartment", "Stage", 1.2, kingdom) +
    theme(legend.position = "none")

op3 <- pcoa(dmr, des, 34, "Management", "Stage", 1.2, kingdom) +
    theme(legend.position = "none")

###############################################################################
###############################################################################
###############################################################################
## Output the files
library(gridExtra)

all_p123 <- grid.arrange(p1, fp1, op1,
                        p2, fp2, op2,
                        p3, fp3, op3,
                        nrow = 3, ncol = 3)
ggsave("Figure_2a.pdf", all_p123, width = 10, height = 10)
