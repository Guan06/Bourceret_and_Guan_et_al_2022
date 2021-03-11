#!/biodata/dep_psl/grp_rgo/tools/bin/Rscript

library(parallelDist)
source("../00.common_scripts/plot_settings.R")
source("./figure_s10_settings.R")

design_file <- "../00.data/design_48.txt"
design <- read.table(design_file, header = T, sep = "\t")
design$Man_Sta <- paste0(design$Management, "_", design$Stage)

###############################################################################
###############################################################################
## Constrained Bacteria and Fungi PCoA plots
library(scales)
library(vegan)
source("./cpcoa_functions.R")

###############################################################################
## for Bacteria
bac_file <- "../00.data/final_ASV/Bacteria/ASV_raref.rds"
bac <- readRDS(bac_file)
response <- bac[, colnames(bac) %in% design$Sample_ID]
response <- response[rowSums(response) > 0, ]

this_design <- design[design$Sample_ID %in% colnames(response), ]
this_design <- this_design[match(colnames(response), this_design$Sample_ID), ]

response <- as.matrix(t(response))
bc <- as.matrix(parDist(response, method = "bray"))

capscale.gen <- capscale(bc ~ Genotype * Management * Stage,
                         data = this_design, add = F, sqrt.dist = T)
ad <- adonis(formula = bc ~ Management * Stage * Genotype,
             data = this_design, add = F, by = "margin")
print(ad)

# ANOVA-like permutation analysis
perm_anova.gen <- anova.cca(capscale.gen)
p.val <- perm_anova.gen[1, 4]
# generate variability tables and calculate confidence intervals for the
# variance
var_tbl.gen <- variability_table(capscale.gen)
eig <- capscale.gen$CCA$eig

variance <- capscale.gen$CCA$tot.chi / capscale.gen$tot.chi
variance <- variance * 100

if (ncol(capscale.gen$CCA$wa) >= 2) {

    points <- data.frame(x = capscale.gen$CCA$wa[, 1],
                         y = capscale.gen$CCA$wa[, 2])

    this_design <- this_design[match(rownames(points), this_design$Sample_ID), ]
    points <- cbind(points, this_design$Genotype, this_design$Man_Sta)
    colnames(points)[3:4] <- c("Genotype", "Man_Sta")
    points$Man_Sta <- factor(points$Man_Sta, levels = order0, ordered= TRUE)

    eig1 <- 100 * eig[1] / sum(eig)
    eig2 <- 100 * eig[2] / sum(eig)

    p1 <- ggplot(points, aes(x, y, shape = Man_Sta, color = Genotype)) +
                geom_point(size = 3, alpha = 0.8) +
                scale_colour_manual(values = as.character(c_Gen$color)) +
                scale_shape_manual(values = s_Man_Sta$shape) +
                main_theme +
                labs(x = paste0("CPCo 1 (", format(eig1, digits = 4), "%)"),
                     y = paste0("CPCo 2 (", format(eig2, digits = 4), "%)")) +
                ggtitle(paste0("Bac: ",
                               format(variance, digits = 3),
                               " % of variance; p = ",
                               format(p.val, digits = 2))) +
                theme(legend.position = "none",
                      plot.title = element_text(size = 12))
}

###############################################################################
## for Fungi
fun_file <- "../00.data/final_ASV/Fungi/ASV_raref.rds"
fun <- readRDS(fun_file)
response <- fun[, colnames(fun) %in% design$Sample_ID]
response <- response[rowSums(response) > 0, ]

this_design <- design[design$Sample_ID %in% colnames(response), ]
this_design <- this_design[match(colnames(response), this_design$Sample_ID), ]

response <- as.matrix(t(response))
bc <- as.matrix(parDist(response, method = "bray"))

capscale.gen <- capscale(bc ~ Genotype * Management * Stage,
                         data = this_design, add = F, sqrt.dist = T)
ad <- adonis(formula = bc ~ Stage * Management * Genotype,
             data = this_design, add = F, by = "margin")
print(ad)

# ANOVA-like permutation analysis
perm_anova.gen <- anova.cca(capscale.gen)
p.val <- perm_anova.gen[1, 4]
var_tbl.gen <- variability_table(capscale.gen)
eig <- capscale.gen$CCA$eig

variance <- capscale.gen$CCA$tot.chi / capscale.gen$tot.chi
variance <- variance * 100

if (ncol(capscale.gen$CCA$wa) >= 2) {

    points <- data.frame(x = capscale.gen$CCA$wa[, 1],
                         y = capscale.gen$CCA$wa[, 2])

    this_design <- this_design[match(rownames(points), this_design$Sample_ID), ]
    points <- cbind(points, this_design$Genotype, this_design$Man_Sta)
    colnames(points)[3:4] <- c("Genotype", "Man_Sta")
    points$Man_Sta <- factor(points$Man_Sta, levels = order0, ordered= TRUE)

    eig1 <- 100 * eig[1] / sum(eig)
    eig2 <- 100 * eig[2] / sum(eig)

    p2 <- ggplot(points, aes(x, y, shape = Man_Sta, color = Genotype)) +
                geom_point(size = 3, alpha = 0.8) +
                scale_colour_manual(values = as.character(c_Gen$color)) +
                scale_shape_manual(values = s_Man_Sta$shape) +
                main_theme +
                labs(x = paste0("CPCo 1 (", format(eig1, digits = 4), "%)"),
                     y = paste0("CPCo 2 (", format(eig2, digits = 4), "%)")) +
                ggtitle(paste0("Fun: ",
                               format(variance, digits = 3),
                               " % of variance; p = ",
                               format(p.val, digits = 2))) +
                theme(legend.position = "none",
                      plot.title = element_text(size = 12))
}

## put together
library(gridExtra)
all <- grid.arrange(p1, p2, nrow = 1, ncol = 2)
ggsave("Figure_S10b.pdf", all, width = 7, height = 3)
