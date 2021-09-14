#!/biodata/dep_psl/grp_rgo/tools/bin/Rscript

library(ggplot2)
library(ggpubr)

source("../00.common_scripts/plot_settings.R")
source("../00.common_scripts/box_plotting.R")
source("../00.common_scripts/violin_plotting.R")

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

p_c1 <- violin(all, "Compartment", "Stage", "Compartment_Stage", "Shannon",
            size = 2)
p_c1 <- p_c1 + theme(legend.position = "none",
                     axis.text.x = element_blank(),
                     axis.text.y = element_text(size = 18),
                     axis.title.y = element_text(size = 18))

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

write.table(sig12, "Figure_1c_sig_Bacteria.txt",
            quote = F, sep = "\t", row.names = F)

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
p_c2 <- violin(all, "Compartment", "Stage", "Compartment_Stage", "Shannon",
            size = 2)
p_c2<- p_c2 + theme(legend.position = "none",
                    axis.text.x = element_blank(),
                    axis.text.y = element_text(size = 18),
                    axis.title.y = element_text(size = 18))

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

write.table(sig12, "Figure_1c_sig_Fungi.txt", quote = F, sep = "\t", row.names = F)

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

p_c3 <- violin(all, "Compartment", "Stage", "Compartment_Stage", "Shannon",
            size = 2)
p_c3 <- p_c3 + theme(legend.position = "none",
                     axis.text.x = element_blank(),
                     axis.text.y = element_text(size = 18),
                     axis.title.y = element_text(size = 18))

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

write.table(sig12, "Figure_1c_sig_Oomycetes.txt", quote = F, sep = "\t", row.names = F)

