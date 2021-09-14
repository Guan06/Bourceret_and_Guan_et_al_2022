#!/biodata/dep_psl/grp_rgo/tools/bin/Rscript

library(parallelDist)
library(scales)
library(vegan)
library(ggplot2)
source("../00.common_scripts/plot_settings.R")
source("../00.common_scripts/cpcoa_functions.R")
source("./figure_05_settings.R")

design_file <- "../00.data/design_48.txt"
design <- read.table(design_file, header = T, sep = "\t")
design$Man_Sta <- paste0(design$Management, "_", design$Stage)

###############################################################################
## function for meta cPCoA plotting

meta_plot <- function(x, design) {
    response <- x

    rownames(response) <- response$Sample_ID
    response <- response[, -1]
    response <- response[rownames(response) %in% design$Sample_ID, ]

    this_design <- design[design$Sample_ID %in% rownames(response), ]
    this_design <- this_design[match(rownames(response), this_design$Sample_ID), ]

    response <- as.matrix(response)
    bc <- as.matrix(parDist(response, method = "euclidean"))

    capscale.gen <- capscale(bc ~ Genotype * Management * Stage,
                             data = this_design, add = F, sqrt.dist = T)

    ad <- adonis(formula = bc ~  Stage * Genotype * Management,
            data = this_design, add = F, by = "margin")
    print(ad)

    perm_anova.gen <- anova.cca(capscale.gen)
    p.val <- perm_anova.gen[1, 4]
    var_tbl.gen <- variability_table(capscale.gen)
    eig <- capscale.gen$CCA$eig

    variance <- capscale.gen$CCA$tot.chi / capscale.gen$tot.chi
    variance <- variance * 100

    if (ncol(capscale.gen$CCA$wa) < 2) break
    points <- data.frame(x = capscale.gen$CCA$wa[, 1],
                         y = capscale.gen$CCA$wa[, 2])

    this_design <- this_design[match(rownames(points), this_design$Sample_ID), ]
    points <- cbind(points, this_design$Genotype, this_design$Man_Sta)
    colnames(points)[3:4] <- c("Genotype", "Man_Sta")
    points$Man_Sta <- factor(points$Man_Sta, levels = order0, ordered= TRUE)

    eig1 <- 100 * eig[1] / sum(eig)
    eig2 <- 100 * eig[2] / sum(eig)

    p <- ggplot(points, aes(x, y, shape = Man_Sta, color = Genotype)) +
                geom_point(size = 3, alpha = 0.8) +
                scale_colour_manual(values = as.character(c_Gen$color)) +
                scale_shape_manual(values = s_Man_Sta$shape) +
                main_theme +
                labs(x = paste0("CPCo 1 (", format(eig1, digits = 4), "%)"),
                     y = paste0("CPCo 2 (", format(eig2, digits = 4), "%)")) +
                ggtitle(paste0(format(variance, digits = 3),
                               " % of variance; p = ",
                               format(p.val, digits = 2))) +
                theme(legend.position = "none",
                      plot.title = element_text(size = 12))
    return(p)
}

lipid_file <- "../00.data/meta_data/lipid_NK_NPK.txt"
lipid <- read.table(lipid_file, header = T, sep = "\t")
p2 <- meta_plot(lipid, design)
ggsave("Figure_05_a.pdf", p2, width = 3.5, height = 3)
