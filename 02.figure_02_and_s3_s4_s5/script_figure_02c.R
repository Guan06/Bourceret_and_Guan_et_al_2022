#!/biodata/dep_psl/grp_rgo/tools/bin/Rscript

library(ggplot2)
library(ggpubr)

source("../00.common_scripts/plot_settings.R")
source("../00.common_scripts/box_plotting.R")

###############################################################################
###############################################################################
## For Bacteria
alpha_file <- "../00.data/alpha_diversity/alpha_Bacteria.txt"
design_file <- "../00.data/design.txt"
kingdom <- "Bacteria"

alpha <- read.table(alpha_file, header = T, sep = "\t")
design <- read.table(design_file, header = T, sep = "\t")
m <- merge(alpha, design)

###############################################################################
# Grouped by Compartment_Stage, colored by Compartment and shaped by Stage

all <- m[, colnames(m) %in% c("Sample_ID", "Shannon", "Compartment",
                                "Stage", "Compartment_Stage")]

all$Compartment_Stage <- factor(all$Compartment_Stage,
                                levels = Compartment_Stage_levels,
                                ordered = TRUE)

p1 <- box(all, "Compartment", "Stage", "Compartment_Stage", "Shannon")
p1 <- p1 + ggtitle(kingdom) + theme(legend.position = "none")

sig1 <- box_sig(all, "Compartment", "Shannon")
Com_lst <- as.character(unique(all$Compartment))
sig2 <- c()
for (c in Com_lst) {
    this_com <- all[all$Compartment == c, ]
    this_c_sig <- box_sig(this_com, "Stage", "Shannon")
    this_c_sig$Group1 <- paste0(c, "_", this_c_sig$Group1)
    this_c_sig$Group2 <- paste0(c, "_", this_c_sig$Group2)
    sig2 <- rbind(sig2, this_c_sig)
}
sig12 <- rbind(sig1, sig2)
sig12$Significance <- as.numeric(as.character(sig12$Significance))
sig12$Sig <- ifelse(sig12$Significance < 0.05, TRUE, FALSE)
sig12$FDR <- p.adjust(sig12$Significance, method = "fdr")
sig12$Sig_FDR <- ifelse(as.numeric(as.character(sig12$FDR)) < 0.05,
                        TRUE, FALSE)

write.table(sig12, "sig_Bacteria.txt", quote = F, sep = "\t", row.names = F)

###############################################################################
###############################################################################
## For Fungi

alpha_file <- "../00.data/alpha_diversity/alpha_Fungi.txt"
kingdom <- "Fungi"

alpha <- read.table(alpha_file, header = T, sep = "\t")
design <- read.table(design_file, header = T, sep = "\t")
m <- merge(alpha, design)

###############################################################################
# Grouped by Compartment_Stage, colored by Compartment and shaped by Stage
all <- m[, colnames(m) %in% c("Sample_ID", "Shannon", "Compartment",
                                "Stage", "Compartment_Stage")]

all$Compartment_Stage <- factor(all$Compartment_Stage,
                                levels = Compartment_Stage_levels,
                                ordered = TRUE)
fp1 <- box(all, "Compartment", "Stage", "Compartment_Stage", "Shannon")
fp1 <- fp1 + ggtitle(kingdom) + theme(legend.position = "none")

sig1 <- box_sig(all, "Compartment", "Shannon")
Com_lst <- as.character(unique(all$Compartment))
sig2 <- c()
for (c in Com_lst) {
    this_com <- all[all$Compartment == c, ]
    this_c_sig <- box_sig(this_com, "Stage", "Shannon")
    this_c_sig$Group1 <- paste0(c, "_", this_c_sig$Group1)
    this_c_sig$Group2 <- paste0(c, "_", this_c_sig$Group2)
    sig2 <- rbind(sig2, this_c_sig)
}
sig12 <- rbind(sig1, sig2)
sig12$Significance <- as.numeric(as.character(sig12$Significance))
sig12$Sig <- ifelse(sig12$Significance < 0.05, TRUE, FALSE)
sig12$FDR <- p.adjust(sig12$Significance, method = "fdr")
sig12$Sig_FDR <- ifelse(as.numeric(as.character(sig12$FDR)) < 0.05,
                        TRUE, FALSE)

write.table(sig12, "sig_Fungi.txt", quote = F, sep = "\t", row.names = F)

###############################################################################
###############################################################################
## For Oomycetes

alpha_file <- "../00.data/alpha_diversity/alpha_Oomycetes.txt"
kingdom <- "Oomycetes"

alpha <- read.table(alpha_file, header = T, sep = "\t")
design <- read.table(design_file, header = T, sep = "\t")
m <- merge(alpha, design)

###############################################################################
# Grouped by Compartment_Stage, colored by Compartment and shaped by Stage
all <- m[, colnames(m) %in% c("Sample_ID", "Shannon", "Compartment",
                                "Stage", "Compartment_Stage")]

all$Compartment_Stage <- factor(all$Compartment_Stage,
                                levels = Compartment_Stage_levels,
                                ordered = TRUE)

op1 <- box(all, "Compartment", "Stage", "Compartment_Stage", "Shannon")
op1 <- op1 + ggtitle(kingdom) + theme(legend.position = "none")

sig1 <- box_sig(all, "Compartment", "Shannon")
Com_lst <- as.character(unique(all$Compartment))
sig2 <- c()
for (c in Com_lst) {
    this_com <- all[all$Compartment == c, ]
    this_c_sig <- box_sig(this_com, "Stage", "Shannon")
    this_c_sig$Group1 <- paste0(c, "_", this_c_sig$Group1)
    this_c_sig$Group2 <- paste0(c, "_", this_c_sig$Group2)
    sig2 <- rbind(sig2, this_c_sig)
}
sig12 <- rbind(sig1, sig2)
sig12$Significance <- as.numeric(as.character(sig12$Significance))
sig12$Sig <- ifelse(sig12$Significance < 0.05, TRUE, FALSE)
sig12$FDR <- p.adjust(sig12$Significance, method = "fdr")
sig12$Sig_FDR <- ifelse(as.numeric(as.character(sig12$FDR)) < 0.05,
                        TRUE, FALSE)

write.table(sig12, "sig_Oomycetes.txt", quote = F, sep = "\t", row.names = F)

###############################################################################
###############################################################################
## Output
library(gridExtra)

all_p1 <- grid.arrange(p1, fp1, op1, nrow = 1, ncol = 3)
ggsave("Figure_2c.pdf", all_p1, width = 10, height = 5)
