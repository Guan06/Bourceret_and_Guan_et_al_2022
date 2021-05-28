#!/biodata/dep_psl/grp_rgo/tools/bin/Rscript

library(ggplot2)
library(parallelDist)
library(cowplot)
source("../00.common_scripts/plot_settings.R")
source("../00.common_scripts/pcoa_plotting.R")

###############################################################################
###############################################################################
## For Bacteria

asv_file <- "../00.data/final_ASV/Bacteria/ASV_raref.rds"
design_file <- "../00.data/design.txt"
kingdom <- "Bacteria"

###############################################################################
## For soil samples only
asv <- readRDS(asv_file)
design <- read.table(design_file, header = T, sep = "\t")
design <- design[design$Compartment == "Bulksoil", ]
design$Plot <- as.character(design$Plot)

# Get the intersection of asv and design file
inter <- intersect(as.character(colnames(asv)), as.character(design$Sample_ID))
asv <- asv[, colnames(asv) %in% inter]
design <- design[design$Sample_ID %in% inter, ]

# Bray-Curtis distance
bc <- as.matrix(parDist(t(asv), method = "bray"))
des <- design
dis <- bc
dmr <- cmdscale(dis, k = 4, eig = T)
n <- nrow(dis)
print(n)

## Management + Stage
p_a1 <- pcoa(dmr, des, 12, "Plot", "Stage", 3.4, kingdom) +
    ggtitle('Bacteria') +
    theme(legend.position = "none",
          axis.text = element_text(size = 18),
          axis.title.x = element_text(size = 18),
          axis.title.y = element_text(size = 18))

###############################################################################
## For all compartment samples

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
n <- nrow(dis)

## Compartment + Stage / Management + Stage
p_b1 <- pcoa(dmr, des, 12, "Compartment", "Stage", 2.4, kingdom) +
    theme(legend.position = "none",
          axis.text = element_text(size = 18),
          axis.title.x = element_text(size = 18),
          axis.title.y = element_text(size = 18))

p_s4_a1 <- pcoa(dmr, des, 13, "Compartment", "Stage", 2.4, kingdom) +
    ggtitle('Bacteria') +
    theme(legend.position = "none",
          axis.text = element_text(size = 18),
          axis.title.x = element_text(size = 18),
          axis.title.y = element_text(size = 18))

p_s4_b1 <- pcoa(dmr, des, 34, "Management", "Stage", 2.4, kingdom) +
    theme(legend.position = "none",
          axis.text = element_text(size = 18),
          axis.title.x = element_text(size = 18),
          axis.title.y = element_text(size = 18))

##############################################################################
###############################################################################
## For Fungi

asv_file <- "../00.data/final_ASV/Fungi/ASV_raref.rds"
kingdom <- "Fungi"
###############################################################################
## For soil samples only
asv <- readRDS(asv_file)
design <- read.table(design_file, header = T, sep = "\t")
design$Plot <- as.character(design$Plot)
design <- design[design$Compartment == "Bulksoil", ]

# Get the intersection of asv and design file
inter <- intersect(as.character(colnames(asv)), as.character(design$Sample_ID))
asv <- asv[, colnames(asv) %in% inter]
design <- design[design$Sample_ID %in% inter, ]

# Bray-Curtis distance
bc <- as.matrix(parDist(t(asv), method = "bray"))
des <- design
dis <- bc
dmr <- cmdscale(dis, k = 4, eig = T)
n <- nrow(dis)
print(n)

## Management + Stage
p_a2 <- pcoa(dmr, des, 12, "Plot", "Stage", 3.4, kingdom) +
    ggtitle('Fungi') +
    theme(legend.position = "none",
          axis.text = element_text(size = 18),
          axis.title.x = element_text(size = 18),
          axis.title.y = element_text(size = 18))

###############################################################################
## For all compartment
asv_file <- "../00.data/final_ASV/Fungi/ASV_raref.rds"
kingdom <- "Fungi"

asv <- readRDS(asv_file)
design <- read.table(design_file, header = T, sep = "\t")

inter <- intersect(as.character(colnames(asv)), as.character(design$Sample_ID))
asv <- asv[, colnames(asv) %in% inter]
design <- design[design$Sample_ID %in% inter, ]

# Bray-Curtis distance
bc <- as.matrix(parDist(t(asv), method = "bray"))

des <- design
dis <- bc
dmr <- cmdscale(dis, k = 4, eig = T)
n <- nrow(dis)
print(n)

## Comfpartment + Stage / Management + Stage
p_b2 <- pcoa(dmr, des, 12, "Compartment", "Stage", 2.4, kingdom) +
    theme(legend.position = "none",
          axis.text = element_text(size = 18),
          axis.title.x = element_text(size = 18),
          axis.title.y = element_text(size = 18))

p_s4_a2 <- pcoa(dmr, des, 13, "Compartment", "Stage", 2.4, kingdom) +
    ggtitle('Fungi') +
    theme(legend.position = "none",
          axis.text = element_text(size = 18),
          axis.title.x = element_text(size = 18),
          axis.title.y = element_text(size = 18))

p_s4_b2 <- pcoa(dmr, des, 34, "Management", "Stage", 2.4, kingdom) +
    theme(legend.position = "none",
          axis.text = element_text(size = 18),
          axis.title.x = element_text(size = 18),
          axis.title.y = element_text(size = 18))

###############################################################################
###############################################################################
## For Oomycetes

asv_file <- "../00.data/final_ASV/Oomycetes/ASV_raref.rds"
kingdom <- "Oomycetes"

###############################################################################
## For soil samples only
asv <- readRDS(asv_file)
design <- read.table(design_file, header = T, sep = "\t")
design$Plot <- as.character(design$Plot)
design <- design[design$Compartment == "Bulksoil", ]

inter <- intersect(as.character(colnames(asv)), as.character(design$Sample_ID))
asv <- asv[, colnames(asv) %in% inter]
design <- design[design$Sample_ID %in% inter, ]

# Bray-Curtis distance
bc <- as.matrix(parDist(t(asv), method = "bray"))
des <- design
dis <- bc
dmr <- cmdscale(dis, k = 4, eig = T)
n <- nrow(dis)
print(n)

## Management + Stage
p_a3 <- pcoa(dmr, des, 12, "Plot", "Stage", 3.4, kingdom) +
    ggtitle('Oomycetes') +
    theme(legend.position = "none",
          axis.text = element_text(size = 18),
          axis.title.x = element_text(size = 18),
          axis.title.y = element_text(size = 18))

###############################################################################
## For all samples
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
n <- nrow(dis)
print(n)

## Compartment + Stage / Management + Stage
p_b3 <- pcoa(dmr, des, 12, "Management", "Stage", 2.4, kingdom) +
    theme(legend.position = "none",
          axis.text = element_text(size = 18),
          axis.title.x = element_text(size = 18),
          axis.title.y = element_text(size = 18))

p_s4_a3 <- pcoa(dmr, des, 13, "Compartment", "Stage", 2.4, kingdom) +
    ggtitle('Oomycetes') +
    theme(legend.position = "none",
          axis.text = element_text(size = 18),
          axis.title.x = element_text(size = 18),
          axis.title.y = element_text(size = 18))

p_s4_b3 <- pcoa(dmr, des, 34, "Compartment", "Stage", 2.4, kingdom) +
    theme(legend.position = "none",
          axis.text = element_text(size = 18),
          axis.title.x = element_text(size = 18),
          axis.title.y = element_text(size = 18))

source("./script_figure_01c.R")

###############################################################################
## for Figure S4 (c)
adonis <- read.table("../00.data/adonis.txt", header = T, sep = "\t")

## For Bacteria
adonis$Factor <- factor(adonis$Factor, levels = adonis$Factor)

p_s4_c1 <- ggplot(adonis, aes(x = Factor, y = Bacteria)) + main_theme +
    geom_bar(stat = "identity", fill = "gray61") +
    theme(axis.text.y = element_text(size = 18),
          axis.text.x = element_text(angle = 90, size = 6,
                                     vjust = 0.5, hjust = 1),
          axis.title.y = element_text(size = 18)) +
    labs(x = "", y = "") 
## For Fungi
p_s4_c2 <- ggplot(adonis, aes(x = Factor, y = Fungi)) + main_theme +
    geom_bar(stat = "identity", fill = "gray61") +
    theme(axis.text.y = element_text(size = 18),
          axis.text.x = element_text(angle = 90, size = 6,
                                     vjust = 0.5, hjust = 1),
          axis.title.y = element_text(size = 18)) +
    labs(x = "", y = "")
## For Oomycetes
p_s4_c3 <- ggplot(adonis, aes(x = Factor, y = Oomycetes)) + main_theme +
    geom_bar(stat = "identity", fill = "gray61") +
    theme(axis.text.y = element_text(size = 18),
          axis.text.x = element_text(angle = 90, size = 6,
                                     vjust = 0.5, hjust = 1),
          axis.title.y = element_text(size = 18)) +
    labs(x = "", y = "")
###############################################################################
fig_1 <- plot_grid(p_a1, p_a2, p_a3,
                   p_b1, p_b2, p_b3,
                   p_c1, p_c2, p_c3,
                   ncol = 3, nrow = 3,
                   align = "v", axis = "l",
                   rel_heights = c(1, 0.95, 0.75))

ggsave("Figure_1.pdf", fig_1, width = 15, height = 12.5)

fig_s4 <- plot_grid(p_s4_a1, p_s4_a2, p_s4_a3,
                    p_s4_b1, p_s4_b2, p_s4_b3,
                    p_s4_c1, p_s4_c2, p_s4_c3,
                    align = "v", axis = "l",
                    ncol = 3, nrow = 3, rel_heights = c(1, 0.95, 0.75))
ggsave("Figure_S4.pdf", fig_s4, width = 15, height = 12.5)
