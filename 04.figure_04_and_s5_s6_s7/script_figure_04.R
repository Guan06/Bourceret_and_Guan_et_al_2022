#!/biodata/dep_psl/grp_rgo/guan/bin/Rscript

library(parallelDist)
library(vegan)
library(cowplot)
source("../00.common_scripts/plot_settings.R")
source("cpcoa_functions.R")

design_file <- "../00.data/design.txt"
meta_file1 <- "../00.data/meta_data/metabolites_lipid.txt"
meta_file2 <- "../00.data/meta_data/metabolites_aa.txt"
meta_file3 <- "../00.data/meta_data/ionomic_20.txt"

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
                scale_colour_manual(values =  c_Man) +
                scale_shape_manual(values = s_Sta) +
                main_theme +
#                coord_fixed(ratio = 1) +
                labs(x = paste0("CPCo 1 (", format(eig1, digits = 4), "%)"),
                     y = paste0("CPCo 2 (", format(eig2, digits = 4), "%)")) +
                ggtitle(paste0(format(variance, digits = 4),
                               " %; p = ",
                               format(p.val, digits = 2))) +
                theme(legend.position = "none",
                      plot.title = element_text(size = 14))
}

## for lipid
meta <- read.table(meta_file1, header = T, sep="\t")
rownames(meta) <- meta$Sample_ID
meta <- meta[, -1]
meta[is.na(meta)] <- 0

## for each genotype
design <- read.table(design_file, header = T, sep = "\t")
design <- design[design$Genotype != "5_pht1;6", ]
design <- design[design$Sample_ID %in% rownames(meta), ]

g_a1 <- plot_g("1_B73")
g_a2 <- plot_g("3_PH207")
g_a3 <- plot_g("2_DK105")
g_a4 <- plot_g("4_F2")

## for aa
meta <- read.table(meta_file2, header = T, sep="\t")
rownames(meta) <- meta$Sample_ID
meta <- meta[, -1]
meta[is.na(meta)] <- 0

design <- read.table(design_file, header = T, sep = "\t")
design <- design[design$Genotype != "5_pht1;6", ]
design <- design[design$Sample_ID %in% rownames(meta), ]

g_b1 <- plot_g("1_B73")
g_b2 <- plot_g("3_PH207")
g_b3 <- plot_g("2_DK105")
g_b4 <- plot_g("4_F2")

## for ion
meta <- read.table(meta_file3, header = T, sep="\t")
rownames(meta) <- meta$Sample_ID
meta <- meta[, -1]
meta[is.na(meta)] <- 0

design <- read.table(design_file, header = T, sep = "\t")
design <- design[design$Genotype != "5_pht1;6", ]
design <- design[design$Sample_ID %in% rownames(meta), ]

g_c1 <- plot_g("1_B73")
g_c2 <- plot_g("3_PH207")
g_c3 <- plot_g("2_DK105")
g_c4 <- plot_g("4_F2")

all <- plot_grid(g_a1, g_a2, g_a3, g_a4,
                 g_b1, g_b2, g_b3, g_b4,
                 g_c1, g_c2, g_c3, g_c4,
                 align = "vh", ncol = 4, nrow = 3)
ggsave("Figure_4.pdf", all, width = 12, height = 9)
