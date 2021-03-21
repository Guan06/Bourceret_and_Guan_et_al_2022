#!/biodata/dep_psl/grp_rgo/guan/bin/Rscript

library(parallelDist)
library(vegan)
library(gridExtra)
source("../00.common_scripts/plot_settings.R")
source("cpcoa_functions.R")

args = commandArgs(trailingOnly = TRUE)
design_file <- args[1]
meta_file <- args[2]
plot_file <- args[3]

design <- read.table(design_file, header = T, sep = "\t")
design <- design[design$Genotype != "5_pht1;6", ]

meta <- read.table(meta_file, header = T, sep="\t")
rownames(meta) <- meta$Sample_ID
meta <- meta[, -1]
meta[is.na(meta)] <- 0

## for each genotype
design <- design[design$Sample_ID %in% rownames(meta), ]

plot_g <- function(g) {
    this_des <- design[design$Genotype == g, ]
    this_meta <- meta[rownames(meta) %in% this_des$Sample_ID, ]
    this_meta_scale <- apply(this_meta, 2, function(x)
                       scale(x, center = TRUE, scale = TRUE))
    rownames(this_meta_scale) <- rownames(this_meta)
    # Euclidean distance between samples
    dis <- as.matrix(parDist(as.matrix(this_meta_scale)))
    
    capscale.gen <- capscale(dis ~ Management * Stage,
                             data = this_des, add = F, sqrt.dist = T)
    perm_anova.gen <- anova.cca(capscale.gen)
    p.val <- perm_anova.gen[1, 4]
    var_tbl.gen <- variability_table(capscale.gen)
    eig <- capscale.gen$CCA$eig

    variance <- capscale.gen$CCA$tot.chi / capscale.gen$tot.chi
    variance <- variance * 100

    if (ncol(capscale.gen$CCA$wa) < 2) break
    points <- data.frame(x = capscale.gen$CCA$wa[, 1],
                         y = capscale.gen$CCA$wa[, 2])

    this_des <- this_des[match(rownames(points), this_des$Sample_ID), ]
    points <- cbind(points, this_des)

    eig1 <- 100 * eig[1] / sum(eig)
    eig2 <- 100 * eig[2] / sum(eig)

    p <- ggplot(points, aes(x, y, shape = Stage, color = Management)) +
                geom_point(size = 3, alpha = 0.8) +
                scale_colour_manual(values = c(c_green, c_red, c_grey, c_black)) +
                scale_shape_manual(values = c(16, 1)) +
                main_theme +
#                coord_fixed(ratio = 1) +
                labs(x = paste0("CPCo 1 (", format(eig1, digits = 4), "%)"),
                     y = paste0("CPCo 2 (", format(eig2, digits = 4), "%)")) +
                ggtitle(paste0(format(variance, digits = 3),
                               " %; p = ",
                               format(p.val, digits = 2))) +
                theme(legend.position = "none")
}

g1 <- plot_g("1_B73")
g2 <- plot_g("3_PH207")
g3 <- plot_g("2_DK105")
g4 <- plot_g("4_F2")

all <- grid.arrange(g1, g2, g3, g4, nrow = 1, ncol = 4)
ggsave(plot_file, all, width = 20, height = 5)
