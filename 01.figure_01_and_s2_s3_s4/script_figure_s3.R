#!/biodata/dep_psl/grp_rgo/tools/bin/Rscript

library(ggplot2)
library(cowplot)
source("../00.common_scripts/plot_settings.R")
source("../00.common_scripts/box_plotting.R")

###############################################################################
asv_file <- "../00.data/final_ASV/Bacteria/ASV_raref.rds"
alpha_file <- "../00.data/alpha_diversity/alpha_Bacteria.txt"
design_file <- "../00.data/design.txt"
kingdom <- "Bacteria"

asv <- readRDS(asv_file)
design <- read.table(design_file, header = T, sep = "\t")
design <- design[design$Compartment == "Bulksoil", ]
design$Plot <- as.character(design$Plot)

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
                                labels = Plot_Stage_labels,
                                ordered = TRUE)

p0 <- box(alpha_soil, "Plot", "Stage", "Plot_Stage", "Shannon", size = 1.8) +
    theme(legend.position = "none",
          axis.text.x = element_text(size = 12),
          axis.text.y = element_text(size = 18)) +
    labs(y = '')

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
write.table(sig12, "Figure_s3_sig_Bacteria.txt", quote = F, sep = "\t", row.names = F)

###############################################################################
## For Fungi

asv_file <- "../00.data/final_ASV/Fungi/ASV_raref.rds"
alpha_file <- "../00.data/alpha_diversity/alpha_Fungi.txt"
kingdom <- "Fungi"

asv <- readRDS(asv_file)
design <- read.table(design_file, header = T, sep = "\t")
design$Plot <- as.character(design$Plot)
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
                                labels = Plot_Stage_labels,
                                ordered = TRUE)

fp0 <- box(alpha_soil, "Plot", "Stage", "Plot_Stage", "Shannon", size = 1.8) +
    theme(legend.position = "none",
          axis.text.x = element_text(size = 12),
          axis.text.y = element_text(size = 18)) +
    labs(y = '')

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

write.table(sig12, "Figure_s3_sig_Fungi.txt", quote = F, sep = "\t", row.names = F)

###############################################################################
## For Oomycetes

asv_file <- "../00.data/final_ASV/Oomycetes/ASV_raref.rds"
alpha_file <- "../00.data/alpha_diversity/alpha_Oomycetes.txt"
kingdom <- "Oomycetes"

asv <- readRDS(asv_file)
design <- read.table(design_file, header = T, sep = "\t")
design$Plot <- as.character(design$Plot)
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
                                labels = Plot_Stage_labels,
                                ordered = TRUE)

op0 <- box(alpha_soil, "Plot", "Stage", "Plot_Stage", "Shannon", size = 1.8) +
    theme(legend.position = "none",
          axis.text.x = element_text(size = 12),
          axis.text.y = element_text(size = 18)) +
    labs(y = '')

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

write.table(sig12, "Figure_s3_sig_Oomycetes.txt", quote = F, sep = "\t", row.names = F)

###############################################################################
## Output the files
all <- plot_grid(p0, fp0, op0, ncol = 3, align = "h", labels = 'auto')
ggsave("Figure_S3.pdf", all, width = 14, height = 4)
