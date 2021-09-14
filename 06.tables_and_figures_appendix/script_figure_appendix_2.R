#!/biodata/dep_psl/grp_rgo/tools/bin/Rscript
library(ggplot2)
library(parallelDist)
library(ggpubr)
source("../00.common_scripts/plot_settings.R")
source("../00.common_scripts/pcoa_plotting.R")
source("../00.common_scripts/box_plotting.R")
source("../00.common_scripts/violin_plotting.R")

order_new <- c("1_B73_Vegetative",
               "1_B73_Reproductive",
               "3_PH207_Vegetative",
               "3_PH207_Reproductive",
               "5_pht1;6_Vegetative",
               "5_pht1;6_Reproductive",
               "2_DK105_Vegetative",
               "2_DK105_Reproductive",
               "4_F2_Vegetative",
               "4_F2_Reproductive")
label_new <- c("B73_VG", "B73_RP", "PH207_VG", "PH207_RP",
               "pht1;6_Vg", "pht1;6_RP",
               "DK105_VG", "DK105_RP", "F2_VG", "F2_RP")

###############################################################################
# Bacteria
asv_file <- "../00.data/final_ASV/Bacteria/ASV_raref.rds"
design_file <- "../00.data/design.txt"
alpha_file <- "../00.data/alpha_diversity/alpha_Bacteria.txt"
kingdom <- "Bacteria"

## For Root samples PCoA plots
asv <- readRDS(asv_file)
design <- read.table(design_file, header = T, sep = "\t")
design <- design[design$Compartment == "Root", ]

inter <- intersect(as.character(colnames(asv)), as.character(design$Sample_ID))
asv <- asv[, colnames(asv) %in% inter]
design <- design[design$Sample_ID %in% inter, ]

dis <- as.matrix(parDist(t(asv), method = "bray"))
dmr <- cmdscale(dis, k = 4, eig = T)

p1a <- pcoa(dmr, design, 12, "Genotype", "Stage", 2.6, kingdom) +
        theme(legend.position = "none", axis.text = element_text(size = 18),
        axis.title.x = element_text(size = 18),
        axis.title.y = element_text(size = 18))

## For root samples box plots
alpha <- read.table(alpha_file, header = T, sep = "\t")
alpha <- alpha[alpha$Sample_ID %in% design$Sample_ID, ]

design$Pool <- ifelse(design$Genotype %in% c("2_DK105", "4_F2"),
                      "Flint", "Dent")

root <- merge(alpha, design)

root$Genotype_Stage <- paste0(root$Genotype, "_", root$Stage)
root$Genotype_Stage <- factor(root$Genotype_Stage,
                              levels = order_new,
                              ordered = TRUE)
p1b <- violin(root, "Genotype", "Stage", "Genotype_Stage", "Shannon") +
   theme(legend.position = "none",
          axis.title.y = element_text(size = 18),
          axis.text = element_text(size = 18)) +
    scale_x_discrete(limits = order_new, labels = label_new)

lst_genotype <- as.character(unique(root$Genotype))
sig1 <- box_sig(root, "Pool", "Shannon")
sig2 <- box_sig(root, "Genotype", "Shannon")
sig0 <- c()
for (mn in lst_genotype) {
    this <- root[root$Genotype == mn, ]
    this_sig <- box_sig(this, "Stage", "Shannon")
    this_sig$Group1 <- paste0(mn, "_", this_sig$Group1)
    this_sig$Group2 <- paste0(mn, "_", this_sig$Group2)
    sig0 <- rbind(sig0, this_sig)
}
sig_bac_root <- rbind(sig1, sig2, sig0)
sig_bac_root$Significance <- as.numeric(as.character(sig_bac_root$Significance))
sig_bac_root$Sig <- ifelse(sig_bac_root$Significance < 0.05, TRUE, FALSE)
sig_bac_root$FDR <- p.adjust(sig_bac_root$Significance, method = "fdr")
sig_bac_root$Sig_FDR <- ifelse(as.numeric(as.character(sig_bac_root$FDR)) < 0.05,
                        TRUE, FALSE)

write.table(sig_bac_root, "Figure_appendix_2_sig_bac_root.txt",
            quote = F, sep = "\t", row.names = F)

###############################################################################
# Fungi
asv_file <- "../00.data/final_ASV/Fungi/ASV_raref.rds"
alpha_file <- "../00.data/alpha_diversity/alpha_Fungi.txt"
kingdom <- "Fungi"

## For Root samples PCoA plots
asv <- readRDS(asv_file)
design <- read.table(design_file, header = T, sep = "\t")
design <- design[design$Compartment == "Root", ]

inter <- intersect(as.character(colnames(asv)), as.character(design$Sample_ID))
asv <- asv[, colnames(asv) %in% inter]
design <- design[design$Sample_ID %in% inter, ]

dis <- as.matrix(parDist(t(asv), method = "bray"))
dmr <- cmdscale(dis, k = 4, eig = T)

p2a <- pcoa(dmr, design, 12, "Genotype", "Stage", 2.6, kingdom) +
    theme(legend.position = "none", axis.text = element_text(size = 18),
        axis.title.x = element_text(size = 18),
        axis.title.y = element_text(size = 18))

## For root samples box plots
alpha <- read.table(alpha_file, header = T, sep = "\t")
alpha <- alpha[alpha$Sample_ID %in% design$Sample_ID, ]

design$Pool <- ifelse(design$Genotype %in% c("2_DK105", "4_F2"),
                      "Flint", "Dent")

root <- merge(alpha, design)

root$Genotype_Stage <- paste0(root$Genotype, "_", root$Stage)
root$Genotype_Stage <- factor(root$Genotype_Stage,
                            levels = order_new,
                            ordered = TRUE)
p2b <- violin(root, "Genotype", "Stage", "Genotype_Stage", "Shannon") +
    theme(legend.position = "none",
          axis.title.y = element_text(size = 18),
          axis.text = element_text(size = 18)) +
    scale_x_discrete(limits = order_new, labels = label_new)

lst_genotype <- as.character(unique(root$Genotype))
sig1 <- box_sig(root, "Pool", "Shannon")
sig2 <- box_sig(root, "Genotype", "Shannon")
sig0 <- c()
for (mn in lst_genotype) {
    this <- root[root$Genotype == mn, ]
    this_sig <- box_sig(this, "Stage", "Shannon")
    this_sig$Group1 <- paste0(mn, "_", this_sig$Group1)
    this_sig$Group2 <- paste0(mn, "_", this_sig$Group2)
    sig0 <- rbind(sig0, this_sig)
}

sig_fun_root <- rbind(sig1, sig2, sig0)
sig_fun_root$Significance <- as.numeric(as.character(sig_fun_root$Significance))
sig_fun_root$Sig <- ifelse(sig_fun_root$Significance < 0.05, TRUE, FALSE)
sig_fun_root$FDR <- p.adjust(sig_fun_root$Significance, method = "fdr")
sig_fun_root$Sig_FDR <- ifelse(as.numeric(as.character(sig_fun_root$FDR)) < 0.05,
                        TRUE, FALSE)

write.table(sig_fun_root, "Figure_appendix_2_sig_fun_root.txt",
            quote = F, sep = "\t", row.names = F)

###############################################################################
# Oomycetes
asv_file <- "../00.data/final_ASV/Oomycetes/ASV_raref.rds"
alpha_file <- "../00.data/alpha_diversity/alpha_Oomycetes.txt"
kingdom <- "Oomycetes"

## For Root samples PCoA plots
asv <- readRDS(asv_file)
design <- read.table(design_file, header = T, sep = "\t")
design <- design[design$Compartment == "Root", ]

inter <- intersect(as.character(colnames(asv)), as.character(design$Sample_ID))
asv <- asv[, colnames(asv) %in% inter]
design <- design[design$Sample_ID %in% inter, ]

dis <- as.matrix(parDist(t(asv), method = "bray"))
dmr <- cmdscale(dis, k = 4, eig = T)

p3a <- pcoa(dmr, design, 12, "Genotype", "Stage", 2, kingdom) +
    theme(legend.position = "none", axis.text = element_text(size = 18),
        axis.title.x = element_text(size = 18),
        axis.title.y = element_text(size = 18))

## For root samples box plots
alpha <- read.table(alpha_file, header = T, sep = "\t")
alpha <- alpha[alpha$Sample_ID %in% design$Sample_ID, ]

design$Pool <- ifelse(design$Genotype %in% c("2_DK105", "4_F2"),
                      "Flint", "Dent")
root <- merge(alpha, design)

root$Genotype_Stage <- paste0(root$Genotype, "_", root$Stage)
root$Genotype_Stage <- factor(root$Genotype_Stage,
                            levels = order_new,
                            ordered = TRUE)

p3b <- violin(root, "Genotype", "Stage", "Genotype_Stage", "Shannon") +
    theme(legend.position = "none",
          axis.title.y = element_text(size = 18),
          axis.text = element_text(size = 18)) +
    scale_x_discrete(limits = order_new, labels = label_new)

lst_genotype <- as.character(unique(root$Genotype))
sig1 <- box_sig(root, "Pool", "Shannon")
sig2 <- box_sig(root, "Genotype", "Shannon")
sig0 <- c()
for (mn in lst_genotype) {
    this <- root[root$Genotype == mn, ]
    this_sig <- box_sig(this, "Stage", "Shannon")
    this_sig$Group1 <- paste0(mn, "_", this_sig$Group1)
    this_sig$Group2 <- paste0(mn, "_", this_sig$Group2)
    sig0 <- rbind(sig0, this_sig)
}
sig_oom_root <- rbind(sig1, sig2, sig0)
sig_oom_root$Significance <- as.numeric(as.character(sig_oom_root$Significance))
sig_oom_root$Sig <- ifelse(sig_oom_root$Significance < 0.05, TRUE, FALSE)
sig_oom_root$FDR <- p.adjust(sig_oom_root$Significance, method = "fdr")
sig_oom_root$Sig_FDR <- ifelse(as.numeric(as.character(sig_oom_root$FDR)) < 0.05,
                        TRUE, FALSE)

write.table(sig_oom_root, "Figure_appendix_2_sig_oom_root.txt",
            quote = F, sep = "\t", row.names = F)

## combine together
library(cowplot)
all <- plot_grid(p1a, p2a, p3a,
        p1b, p2b, p3b,
        nrow = 2, ncol = 3,
        align = "v", axis = "l",rel_heights = c(1, 0.75))

ggsave("Figure_appendix_2.pdf", all, width = 16, height = 8.5)
