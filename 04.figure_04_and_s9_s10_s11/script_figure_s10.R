#!/biodata/dep_psl/grp_rgo/guan/bin/Rscript

library(parallelDist)
library(vegan)
library(gridExtra)
source("../00.common_scripts/plot_settings.R")
source("../00.common_scripts/cpcoa_functions.R")

design_file <- "../00.data/design.txt"
sugar_file <- "../00.data/meta_data/metabolites_sugar.txt"

design <- read.table(design_file, header = T, sep = "\t")
design <- design[design$Field == "DEMO", ]

sugar <- read.table(sugar_file, header = T, sep="\t")
rownames(sugar) <- sugar$Sample_ID
sugar <- sugar[, -1]

design <- design[design$Sample_ID %in% rownames(sugar), ]
glst <- unique(design$Genotype)

plot_g <- function(g){
    this_des <- design[design$Genotype == g, ]
    this_meta <- sugar[rownames(sugar) %in% this_des$Sample_ID, ]
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
                scale_colour_manual(values = c_Man) +
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

g1 <- plot_g("1_B73")
g2 <- plot_g("3_PH207")
g3 <- plot_g("2_DK105")
g4 <- plot_g("4_F2")

all <- grid.arrange(g1, g2, g3, g4, nrow = 1, ncol = 4)
ggsave("Figure_S10.pdf", all, width = 12, height = 3)
