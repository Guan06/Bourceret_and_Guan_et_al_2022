#!/biodata/dep_psl/grp_rgo/tools/bin/Rscript

library(ggplot2)
library(parallelDist)
source("../00.common_scripts/plot_settings.R")
source("../00.common_scripts/pcoa_plotting.R")
source("../00.common_scripts/box_plotting.R")

###############################################################################
###############################################################################
## For Bacteria

asv_file <- "../00.data/final_ASV/Bacteria/ASV_raref.rds"
alpha_file <- "../00.data/alpha_diversity/alpha_Bacteria.txt"
design_file <- "../00.data/design.txt"
kingdom <- "Bacteria"

asv <- readRDS(asv_file)
design <- read.table(design_file, header = T, sep = "\t")
design <- design[design$Compartment == "Bulksoil", ]

alpha <- read.table(alpha_file, header = T, sep = "\t")
alpha_soil <- merge(alpha, design)

# Group alpha diversity of soil samples by Management_Stage
alpha_soil <- alpha_soil[, colnames(alpha_soil) %in% c("Sample_ID", "Shannon",
                                                       "Management", "Plot",
                                                       "Stage", "Field")]
alpha_soil$Management_Stage <- paste0(alpha_soil$Management, "_",
                                      alpha_soil$Stage)

alpha_soil$Management_Stage <- factor(alpha_soil$Management_Stage,
                                      levels = Management_Stage_levels,
                                      ordered = TRUE)

alpha_soil$Plot_Stage <- paste0(alpha_soil$Plot, "_", alpha_soil$Stage)
alpha_soil$Plot_Stage <- factor(alpha_soil$Plot_Stage,
                                levels = Plot_Stage_levels,
                                ordered = TRUE)

p0 <- box(alpha_soil, "Plot", "Stage", "Plot_Stage", "Shannon") +
    theme(legend.position = "none",
          axis.text = element_text(size = 18))

# calculate the significance of alpha diversity
sig0 <- box_sig(alpha_soil, "Field", "Shannon")
sig1 <- box_sig(alpha_soil, "Management", "Shannon")
Plot_lst <- as.character(unique(alpha_soil$Plot))
sig2 <- c()
for (p in Plot_lst) {
    this_alpha_soil <- alpha_soil[alpha_soil$Plot == p, ]
    this_p_sig <- box_sig(this_alpha_soil, "Stage", "Shannon")
    this_p_sig$Group1 <- paste0(p, "_", this_p_sig$Group1)
    this_p_sig$Group2 <- paste0(p, "_", this_p_sig$Group2)
    sig2 <- rbind(sig2, this_p_sig)
}
sig12 <- rbind(sig0, sig1, sig2)
sig12$Significance <- as.numeric(as.character(sig12$Significance))
sig12$Sig <- ifelse(sig12$Significance < 0.05, TRUE, FALSE)
sig12$FDR <- p.adjust(sig12$Significance, method = "fdr")
sig12$Sig_FDR <- ifelse(as.numeric(as.character(sig12$FDR)) < 0.05,
                        TRUE, FALSE)
write.table(sig12, "sig_Bacteria.txt", quote = F, sep = "\t", row.names = F)

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
p1 <- pcoa(dmr, des, 12, "Plot", "Stage", 2.4, kingdom) +
    theme(legend.position = "none",
          axis.text = element_text(size = 18),
          axis.text.x = element_text(size = 18))

##############################################################################
###############################################################################
## For Fungi

asv_file <- "../00.data/final_ASV/Fungi/ASV_raref.rds"
alpha_file <- "../00.data/alpha_diversity/alpha_Fungi.txt"
kingdom <- "Fungi"

asv <- readRDS(asv_file)
design <- read.table(design_file, header = T, sep = "\t")
design <- design[design$Compartment == "Bulksoil", ]

alpha <- read.table(alpha_file, header = T, sep = "\t")
alpha_soil <- merge(alpha, design)

# Group alpha diversity of soil samples by Management_Stage
alpha_soil <- alpha_soil[, colnames(alpha_soil) %in% c("Sample_ID", "Shannon",
                                                       "Management", "Plot",
                                                       "Stage", "Field")]
alpha_soil$Management_Stage <- paste0(alpha_soil$Management, "_",
                                      alpha_soil$Stage)

alpha_soil$Management_Stage <- factor(alpha_soil$Management_Stage,
                                      levels = Management_Stage_levels,
                                      ordered = TRUE)
alpha_soil$Plot_Stage <- paste0(alpha_soil$Plot, "_", alpha_soil$Stage)
alpha_soil$Plot_Stage <- factor(alpha_soil$Plot_Stage,
                                levels = Plot_Stage_levels,
                                ordered = TRUE)

fp0 <- box(alpha_soil, "Plot", "Stage", "Plot_Stage", "Shannon") +
    theme(legend.position = "none",
          axis.text = element_text(size = 18),
          axis.text.x = element_text(size = 18))

# calculate the significance of alpha diversity
sig0 <- box_sig(alpha_soil, "Field", "Shannon")
sig1 <- box_sig(alpha_soil, "Management", "Shannon")
Plot_lst <- as.character(unique(alpha_soil$Plot))
sig2 <- c()
for (p in Plot_lst) {
    this_alpha_soil <- alpha_soil[alpha_soil$Plot == p, ]
    this_p_sig <- box_sig(this_alpha_soil, "Stage", "Shannon")
    this_p_sig$Group1 <- paste0(p, "_", this_p_sig$Group1)
    this_p_sig$Group2 <- paste0(p, "_", this_p_sig$Group2)
    sig2 <- rbind(sig2, this_p_sig)
}
sig12 <- rbind(sig0, sig1, sig2)
sig12$Significance <- as.numeric(as.character(sig12$Significance))
sig12$Sig <- ifelse(sig12$Significance < 0.05, TRUE, FALSE)
sig12$FDR <- p.adjust(sig12$Significance, method = "fdr")
sig12$Sig_FDR <- ifelse(as.numeric(as.character(sig12$FDR)) < 0.05,
                        TRUE, FALSE)

write.table(sig12, "sig_Fungi.txt", quote = F, sep = "\t", row.names = F)

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
## Management + Stage
fp1 <- pcoa(dmr, des, 12, "Plot", "Stage", 2.4, kingdom) +
    theme(legend.position = "none",
          axis.text = element_text(size = 18),
          axis.text.x = element_text(size = 18))

###############################################################################
###############################################################################
## For Oomycetes

asv_file <- "../00.data/final_ASV/Oomycetes/ASV_raref.rds"
alpha_file <- "../00.data/alpha_diversity/alpha_Oomycetes.txt"
kingdom <- "Oomycetes"

asv <- readRDS(asv_file)
design <- read.table(design_file, header = T, sep = "\t")
design <- design[design$Compartment == "Bulksoil", ]

alpha <- read.table(alpha_file, header = T, sep = "\t")
alpha_soil <- merge(alpha, design)

# Group alpha diversity of soil samples by Management_Stage
alpha_soil <- alpha_soil[, colnames(alpha_soil) %in% c("Sample_ID", "Shannon",
                                                       "Management", "Plot",
                                                       "Stage", "Field")]
alpha_soil$Management_Stage <- paste0(alpha_soil$Management, "_",
                                      alpha_soil$Stage)

alpha_soil$Management_Stage <- factor(alpha_soil$Management_Stage,
                                      levels = Management_Stage_levels,
                                      ordered = TRUE)

alpha_soil$Plot_Stage <- paste0(alpha_soil$Plot, "_", alpha_soil$Stage)
alpha_soil$Plot_Stage <- factor(alpha_soil$Plot_Stage,
                                levels = Plot_Stage_levels,
                                ordered = TRUE)

op0 <- box(alpha_soil, "Plot", "Stage", "Plot_Stage", "Shannon") +
    theme(legend.position = "none",
          axis.text = element_text(size = 18),
          axis.text.x = element_text(size = 18))

# calculate the significance of alpha diversity
sig0 <- box_sig(alpha_soil, "Field", "Shannon")
sig1 <- box_sig(alpha_soil, "Management", "Shannon")
Plot_lst <- as.character(unique(alpha_soil$Plot))
sig2 <- c()
for (p in Plot_lst) {
    this_alpha_soil <- alpha_soil[alpha_soil$Plot == p, ]
    this_p_sig <- box_sig(this_alpha_soil, "Stage", "Shannon")
    this_p_sig$Group1 <- paste0(p, "_", this_p_sig$Group1)
    this_p_sig$Group2 <- paste0(p, "_", this_p_sig$Group2)
    sig2 <- rbind(sig2, this_p_sig)
}
sig12 <- rbind(sig0, sig1, sig2)
sig12$Significance <- as.numeric(as.character(sig12$Significance))
sig12$Sig <- ifelse(sig12$Significance < 0.05, TRUE, FALSE)
sig12$FDR <- p.adjust(sig12$Significance, method = "fdr")
sig12$Sig_FDR <- ifelse(as.numeric(as.character(sig12$FDR)) < 0.05,
                        TRUE, FALSE)

write.table(sig12, "sig_Oomycetes.txt", quote = F, sep = "\t", row.names = F)

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
op1 <- pcoa(dmr, des, 12, "Plot", "Stage", 2.4, kingdom) +
    theme(legend.position = "none",
          axis.text = element_text(size = 18),
          axis.text.x = element_text(size = 18))

###############################################################################
###############################################################################
## Output the files
library(gridExtra)

all_p123 <- grid.arrange(p1, fp1, op1,
                         p0, fp0, op0,
                         nrow = 2, ncol = 3)
ggsave("Figure_1.pdf", all_p123, width = 16, height = 10)
