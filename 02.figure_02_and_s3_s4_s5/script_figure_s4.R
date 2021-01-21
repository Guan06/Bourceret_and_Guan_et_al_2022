#!/biodata/dep_psl/grp_rgo/tools/bin/Rscript
library(ggplot2)
library(parallelDist)
library(ggpubr)

source("../00.common_scripts/plot_settings.R")
source("../00.common_scripts/pcoa_plotting.R")
source("../00.common_scripts/violin_plotting.R")
source("../00.common_scripts/box_plotting.R")

order_new <- c("NK_Vegetative",
                "NK_Reproductive",
                "NPK_Vegetative",
                "NPK_Reproductive",
                "CONMIN_Vegetative",
                "CONMIN_Reproductive",
                "BIODYN_Vegetative",
                "BIODYN_Reproductive"
              )

###############################################################################
# Bacteria
## For Rhizosphere samples PCoA plots
asv_file <- "../00.data/final_ASV/Bacteria/ASV_raref.rds"
design_file <- "../00.data/design.txt"
alpha_file <- "../00.data/alpha_diversity/alpha_Bacteria.txt"
kingdom <- "Bacteria"

asv <- readRDS(asv_file)
design <- read.table(design_file, header = T, sep = "\t")
design <- design[design$Compartment == "Rhizosphere", ]

### Get the intersection of asv and design file
inter <- intersect(as.character(colnames(asv)), as.character(design$Sample_ID))
asv <- asv[, colnames(asv) %in% inter]
design <- design[design$Sample_ID %in% inter, ]

### Bray-Curtis distance
dis <- as.matrix(parDist(t(asv), method = "bray"))
dmr <- cmdscale(dis, k = 4, eig = T)
## Compartment + Stage
p1a <- pcoa(dmr, design, 12, "Plot", "Stage", 2, kingdom) +
    theme(legend.position = "none", axis.text = element_text(size = 18)) +
    scale_shape_manual(values = c(1, 16))


alpha <- read.table(alpha_file, header = T, sep = "\t")
alpha <- alpha[alpha$Sample_ID %in% design$Sample_ID, ]
rhi <- merge(alpha, design)

rhi$Management_Stage <- paste0(rhi$Management, "_", rhi$Stage)
rhi$Management_Stage <- factor(rhi$Management_Stage,
                               levels = order_new,
                               ordered = TRUE)

p1b <- violin(rhi, "Management", "Stage", "Management_Stage", "Shannon") +
    theme(legend.position = "none", axis.text = element_text(size = 10)) +
    scale_x_discrete(limits = order_new) +
    scale_shape_manual(values = c(1, 16))


lst_management <- as.character(unique(rhi$Management))
sig_bac_rhi <- box_sig(rhi, "Management", "Shannon")
sig0 <- c()
for (mn in lst_management) {
    this <- rhi[rhi$Management == mn, ]
    this_sig <- box_sig(this, "Stage", "Shannon")
    this_sig$Group1 <- paste0(mn, "_", this_sig$Group1)
    this_sig$Group2 <- paste0(mn, "_", this_sig$Group2)
    sig0 <- rbind(sig0, this_sig)
}
sig_bac_rhi <- rbind(sig0, sig_bac_rhi)
sig_bac_rhi$Significance <- as.numeric(as.character(sig_bac_rhi$Significance))
sig_bac_rhi$Sig <- ifelse(sig_bac_rhi$Significance < 0.05, TRUE, FALSE)
sig_bac_rhi$FDR <- p.adjust(sig_bac_rhi$Significance, method = "fdr")
sig_bac_rhi$Sig_FDR <- ifelse(as.numeric(as.character(sig_bac_rhi$FDR)) < 0.05,
                        TRUE, FALSE)

write.table(sig_bac_rhi, "s4_sig_bac_rhi.txt", quote = F, sep = "\t", row.names = F)

## For Root samples PCoA plots
asv <- readRDS(asv_file)
design <- read.table(design_file, header = T, sep = "\t")
design <- design[design$Compartment == "Root", ]

inter <- intersect(as.character(colnames(asv)), as.character(design$Sample_ID))
asv <- asv[, colnames(asv) %in% inter]
design <- design[design$Sample_ID %in% inter, ]

dis <- as.matrix(parDist(t(asv), method = "bray"))
dmr <- cmdscale(dis, k = 4, eig = T)

p1c <- pcoa(dmr, design, 12, "Plot", "Stage", 2, kingdom) +
        theme(legend.position = "none", axis.text = element_text(size = 18)) +
        scale_shape_manual(values = c(1, 16))

alpha <- read.table(alpha_file, header = T, sep = "\t")
alpha <- alpha[alpha$Sample_ID %in% design$Sample_ID, ]
root <- merge(alpha, design)

root$Management_Stage <- paste0(root$Management, "_", root$Stage)
root$Management_Stage <- factor(root$Management_Stage,
                            levels = order_new,
                            ordered = TRUE)
p1d <- violin(root, "Management", "Stage", "Management_Stage", "Shannon") +
    theme(legend.position = "none", axis.text = element_text(size = 10)) +
    scale_x_discrete(limits = order_new) +
    scale_shape_manual(values = c(1, 16))

lst_management <- as.character(unique(root$Management))
sig_bac_root <- box_sig(root, "Management", "Shannon")
sig0 <- c()
for (mn in lst_management) {
    this <- root[root$Management == mn, ]
    this_sig <- box_sig(this, "Stage", "Shannon")
    this_sig$Group1 <- paste0(mn, "_", this_sig$Group1)
    this_sig$Group2 <- paste0(mn, "_", this_sig$Group2)
    sig0 <- rbind(sig0, this_sig)
}
sig_bac_root <- rbind(sig0, sig_bac_root)
sig_bac_root$Significance <- as.numeric(as.character(sig_bac_root$Significance))
sig_bac_root$Sig <- ifelse(sig_bac_root$Significance < 0.05, TRUE, FALSE)
sig_bac_root$FDR <- p.adjust(sig_bac_root$Significance, method = "fdr")
sig_bac_root$Sig_FDR <- ifelse(as.numeric(as.character(sig_bac_root$FDR)) < 0.05,
                        TRUE, FALSE)

write.table(sig_bac_root, "s4_sig_bac_root.txt", quote = F, sep = "\t", row.names = F)

###############################################################################
# Fungi
asv_file <- "../00.data/final_ASV/Fungi/ASV_raref.rds"
alpha_file <- "../00.data/alpha_diversity/alpha_Fungi.txt"
kingdom <- "Fungi"

## For Rhizosphere samples PCoA plots
asv <- readRDS(asv_file)
design <- read.table(design_file, header = T, sep = "\t")
design <- design[design$Compartment == "Rhizosphere", ]

# Get the intersection of asv and design file
inter <- intersect(as.character(colnames(asv)), as.character(design$Sample_ID))
asv <- asv[, colnames(asv) %in% inter]
design <- design[design$Sample_ID %in% inter, ]

# Bray-Curtis distance
dis <- as.matrix(parDist(t(asv), method = "bray"))
dmr <- cmdscale(dis, k = 4, eig = T)
## Compartment + Stage
p2a <- pcoa(dmr, design, 12, "Plot", "Stage", 2, kingdom) +
    theme(legend.position = "none", axis.text = element_text(size = 18)) +
    scale_shape_manual(values = c(1, 16))


alpha <- read.table(alpha_file, header = T, sep = "\t")
alpha <- alpha[alpha$Sample_ID %in% design$Sample_ID, ]
rhi <- merge(alpha, design)

rhi$Management_Stage <- paste0(rhi$Management, "_", rhi$Stage)
rhi$Management_Stage <- factor(rhi$Management_Stage,
                               levels = order_new,
                               ordered = TRUE)
p2b <- violin(rhi, "Management", "Stage", "Management_Stage", "Shannon") +
    theme(legend.position = "none", axis.text = element_text(size = 10)) +
    scale_x_discrete(limits = order_new) +
    scale_shape_manual(values = c(1, 16))

lst_management <- as.character(unique(rhi$Management))
sig_fun_rhi <- box_sig(rhi, "Management", "Shannon")
sig0 <- c()
for (mn in lst_management) {
    this <- rhi[rhi$Management == mn, ]
    this_sig <- box_sig(this, "Stage", "Shannon")
    this_sig$Group1 <- paste0(mn, "_", this_sig$Group1)
    this_sig$Group2 <- paste0(mn, "_", this_sig$Group2)
    sig0 <- rbind(sig0, this_sig)
}
sig_fun_rhi <- rbind(sig0, sig_fun_rhi)
sig_fun_rhi$Significance <- as.numeric(as.character(sig_fun_rhi$Significance))
sig_fun_rhi$Sig <- ifelse(sig_fun_rhi$Significance < 0.05, TRUE, FALSE)
sig_fun_rhi$FDR <- p.adjust(sig_fun_rhi$Significance, method = "fdr")
sig_fun_rhi$Sig_FDR <- ifelse(as.numeric(as.character(sig_fun_rhi$FDR)) < 0.05,
                        TRUE, FALSE)

write.table(sig_fun_rhi, "s4_sig_fun_rhi.txt", quote = F, sep = "\t", row.names = F)

## For Root samples PCoA plots
asv <- readRDS(asv_file)
design <- read.table(design_file, header = T, sep = "\t")
design <- design[design$Compartment == "Root", ]

inter <- intersect(as.character(colnames(asv)), as.character(design$Sample_ID))
asv <- asv[, colnames(asv) %in% inter]
design <- design[design$Sample_ID %in% inter, ]

dis <- as.matrix(parDist(t(asv), method = "bray"))
dmr <- cmdscale(dis, k = 4, eig = T)

p2c <- pcoa(dmr, design, 12, "Plot", "Stage", 2, kingdom) +
        theme(legend.position = "none", axis.text = element_text(size = 18)) +
        scale_shape_manual(values = c(1, 16))

## For root samples box plots
alpha <- read.table(alpha_file, header = T, sep = "\t")
alpha <- alpha[alpha$Sample_ID %in% design$Sample_ID, ]
root <- merge(alpha, design)

root$Management_Stage <- paste0(root$Management, "_", root$Stage)
root$Management_Stage <- factor(root$Management_Stage,
                            levels = order_new,
                            ordered = TRUE)
p2d <- violin(root, "Management", "Stage", "Management_Stage", "Shannon") +
    theme(legend.position = "none", axis.text = element_text(size = 10)) +
    scale_x_discrete(limits = order_new) +
        scale_shape_manual(values = c(1, 16))

lst_management <- as.character(unique(root$Management))
sig_fun_root <- box_sig(root, "Management", "Shannon")
sig0 <- c()
for (mn in lst_management) {
    this <- root[root$Management == mn, ]
    this_sig <- box_sig(this, "Stage", "Shannon")
    this_sig$Group1 <- paste0(mn, "_", this_sig$Group1)
    this_sig$Group2 <- paste0(mn, "_", this_sig$Group2)
    sig0 <- rbind(sig0, this_sig)
}
sig_fun_root <- rbind(sig0, sig_fun_root)
sig_fun_root$Significance <- as.numeric(as.character(sig_fun_root$Significance))
sig_fun_root$Sig <- ifelse(sig_fun_root$Significance < 0.05, TRUE, FALSE)
sig_fun_root$FDR <- p.adjust(sig_fun_root$Significance, method = "fdr")
sig_fun_root$Sig_FDR <- ifelse(as.numeric(as.character(sig_fun_root$FDR)) < 0.05,
                        TRUE, FALSE)

write.table(sig_fun_root, "s4_sig_fun_root.txt", quote = F, sep = "\t", row.names = F)


###############################################################################
# Oomycetes
asv_file <- "../00.data/final_ASV/Oomycetes/ASV_raref.rds"
alpha_file <- "../00.data/alpha_diversity/alpha_Oomycetes.txt"
kingdom <- "Oomycetes"

## For Rhizosphere samples PCoA plots
asv <- readRDS(asv_file)
design <- read.table(design_file, header = T, sep = "\t")
design <- design[design$Compartment == "Rhizosphere", ]

# Get the intersection of asv and design file
inter <- intersect(as.character(colnames(asv)), as.character(design$Sample_ID))
asv <- asv[, colnames(asv) %in% inter]
design <- design[design$Sample_ID %in% inter, ]

# Bray-Curtis distance
dis <- as.matrix(parDist(t(asv), method = "bray"))
dmr <- cmdscale(dis, k = 4, eig = T)
## Compartment + Stage
p3a <- pcoa(dmr, design, 12, "Plot", "Stage", 2, kingdom) +
    theme(legend.position = "none", axis.text = element_text(size = 18)) +
    scale_shape_manual(values = c(1, 16))

alpha <- read.table(alpha_file, header = T, sep = "\t")
alpha <- alpha[alpha$Sample_ID %in% design$Sample_ID, ]
rhi <- merge(alpha, design)

rhi$Management_Stage <- paste0(rhi$Management, "_", rhi$Stage)
rhi$Management_Stage <- factor(rhi$Management_Stage,
                               levels = order_new,
                               ordered = TRUE)
p3b <- violin(rhi, "Management", "Stage", "Management_Stage", "Shannon") +
    theme(legend.position = "none", axis.text = element_text(size = 10)) +
    scale_x_discrete(limits = order_new) +
    scale_shape_manual(values = c(1, 16))

lst_management <- as.character(unique(rhi$Management))
sig_oom_rhi <- box_sig(rhi, "Management", "Shannon")
sig0 <- c()
for (mn in lst_management) {
    this <- rhi[rhi$Management == mn, ]
    this_sig <- box_sig(this, "Stage", "Shannon")
    this_sig$Group1 <- paste0(mn, "_", this_sig$Group1)
    this_sig$Group2 <- paste0(mn, "_", this_sig$Group2)
    sig0 <- rbind(sig0, this_sig)
}
sig_oom_rhi <- rbind(sig0, sig_oom_rhi)
sig_oom_rhi$Significance <- as.numeric(as.character(sig_oom_rhi$Significance))
sig_oom_rhi$Sig <- ifelse(sig_oom_rhi$Significance < 0.05, TRUE, FALSE)
sig_oom_rhi$FDR <- p.adjust(sig_oom_rhi$Significance, method = "fdr")
sig_oom_rhi$Sig_FDR <- ifelse(as.numeric(as.character(sig_oom_rhi$FDR)) < 0.05,
                        TRUE, FALSE)

write.table(sig_oom_rhi, "s4_sig_oom_rhi.txt", quote = F, sep = "\t", row.names = F)

## For Root samples PCoA plots
asv <- readRDS(asv_file)
design <- read.table(design_file, header = T, sep = "\t")
design <- design[design$Compartment == "Root", ]

inter <- intersect(as.character(colnames(asv)), as.character(design$Sample_ID))
asv <- asv[, colnames(asv) %in% inter]
design <- design[design$Sample_ID %in% inter, ]

dis <- as.matrix(parDist(t(asv), method = "bray"))
dmr <- cmdscale(dis, k = 4, eig = T)

p3c <- pcoa(dmr, design, 12, "Plot", "Stage", 2, kingdom) +
        theme(legend.position = "none", axis.text = element_text(size = 18)) +
        scale_shape_manual(values = c(1, 16))

## For root samples box plots
alpha <- read.table(alpha_file, header = T, sep = "\t")
alpha <- alpha[alpha$Sample_ID %in% design$Sample_ID, ]
root <- merge(alpha, design)

root$Management_Stage <- paste0(root$Management, "_", root$Stage)
root$Management_Stage <- factor(root$Management_Stage,
                            levels = order_new,
                            ordered = TRUE)

p3d <- violin(root, "Management", "Stage", "Management_Stage", "Shannon") +
    theme(legend.position = "none", axis.text = element_text(size = 10)) +
    scale_x_discrete(limits = order_new) +
    scale_shape_manual(values = c(1, 16))

lst_management <- as.character(unique(root$Management))
sig_oom_root <- box_sig(root, "Management", "Shannon")
sig0 <- c()
for (mn in lst_management) {
    this <- root[root$Management == mn, ]
    this_sig <- box_sig(this, "Stage", "Shannon")
    this_sig$Group1 <- paste0(mn, "_", this_sig$Group1)
    this_sig$Group2 <- paste0(mn, "_", this_sig$Group2)
    sig0 <- rbind(sig0, this_sig)
}
sig_oom_root <- rbind(sig0, sig_oom_root)
sig_oom_root$Significance <- as.numeric(as.character(sig_oom_root$Significance))
sig_oom_root$Sig <- ifelse(sig_oom_root$Significance < 0.05, TRUE, FALSE)
sig_oom_root$FDR <- p.adjust(sig_oom_root$Significance, method = "fdr")
sig_oom_root$Sig_FDR <- ifelse(as.numeric(as.character(sig_oom_root$FDR)) < 0.05,
                        TRUE, FALSE)

write.table(sig_oom_root, "s4_sig_oom_root.txt", quote = F, sep = "\t", row.names = F)

## combine together
library(gridExtra)
all <- grid.arrange(p1a, p2a, p3a,
                    p1b, p2b, p3b,
                    p1c, p2c, p3c,
                    p1d, p2d, p3d,
                    nrow = 4, ncol = 3)

ggsave("Figure_S4.pdf", all, width = 16, height = 20)
