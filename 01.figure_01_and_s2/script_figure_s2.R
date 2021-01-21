#!/biodata/dep_psl/grp_rgo/tools/bin/Rscript

library(ggplot2)
library(parallelDist)
source("../00.common_scripts/plot_settings.R")
source("../00.common_scripts/pcoa_plotting.R")
source("../00.common_scripts/box_plotting.R")

design_file <- "../00.data/design.txt"
design <- read.table(design_file, header = T, sep = "\t")
design <- design[design$Stage == "before_sowing", ]
design$Plot <- factor(design$Plot, levels = c_Plo$group)

###############################################################################
## For stable soil properties
soil <- read.table("../00.data/meta_data/soil_properties.txt",
                   header = T, sep="\t")
rownames(soil) <- soil$Sample_ID
soil <- soil[, -1]
soil_scale <- apply(soil, 2, function(x) scale(x, center = TRUE, scale = TRUE))
rownames(soil_scale) <- rownames(soil)

# Euclidean distance between samples based on soil properties
dis <- as.matrix(parDist(as.matrix(soil_scale)))
dmr <- cmdscale(dis, eig = T)
p1 <- pcoa(dmr, design, 12, "Plot", "Stage", 2.4, "Soil properties") +
    theme(legend.position = "none")

## for pH boxplot
ph <- soil[, colnames(soil) %in% c("pH_H2O", "pH_CaCl2")]
ph$pH <- (ph[, 1] + ph[, 2]) / 2
ph$Sample_ID <- rownames(ph)
ph <- merge(ph, design)

ph$Plot <- factor(ph$Plot, levels = c_Plo$group)

p_ph <- ggplot(ph, aes(x = Plot, y = pH)) +
    geom_boxplot(aes(color = Plot), outlier.shape = NA) +
    geom_jitter(aes(color = Plot), size = 2, shape = 2, alpha = 0.8) +
    scale_color_manual(values = as.character(c_Plo$color)) +
    main_theme

sig_ph <- box_sig(ph, "Plot", "pH")
sig_ph$Significance <- as.numeric(as.character(sig_ph$Significance))
sig_ph$Sig <- ifelse(sig_ph$Significance < 0.05, TRUE, FALSE)
sig_ph$FDR <- p.adjust(sig_ph$Significance, method = "fdr")
sig_ph$Sig_FDR <- ifelse(as.numeric(as.character(sig_ph$FDR)) < 0.05,
                        TRUE, FALSE)

write.table(sig_ph, "sig_pH.txt", quote = F, row.names = F, sep = "\t")

## for Clay boxplot
clay <- data.frame(Sample_ID = rownames(soil),
                   Clay = soil$Clay)
clay <- merge(clay, design)
clay$Plot <- factor(clay$Plot, levels = c_Plo$group)

p_clay <- ggplot(clay, aes(x = Plot, y = Clay)) +
    geom_boxplot(aes(color = Plot), outlier.shape = NA) +
    geom_jitter(aes(color = Plot), size = 2, shape = 2, alpha = 0.8) +
    scale_color_manual(values = as.character(c_Plo$color)) +
    main_theme

sig_clay <- box_sig(clay, "Plot", "Clay")
sig_clay$Significance <- as.numeric(as.character(sig_clay$Significance))
sig_clay$Sig <- ifelse(sig_clay$Significance < 0.05, TRUE, FALSE)
sig_clay$FDR <- p.adjust(sig_clay$Significance, method = "fdr")
sig_clay$Sig_FDR <- ifelse(as.numeric(as.character(sig_clay$FDR)) < 0.05,
                        TRUE, FALSE)

write.table(sig_clay, "sig_clay.txt", quote = F, sep = "\t", row.names = F)

###############################################################################
################################# For bacterial community
asv_file <- "../00.data/final_ASV/Bacteria/ASV_raref.rds"
kingdom <- "Bacteria"

asv <- readRDS(asv_file)
asv <- asv[, colnames(asv) %in% design$Sample_ID]

# Bray-Curtis distance
dis <- as.matrix(parDist(t(asv), method = "bray"))
dmr <- cmdscale(dis, k = 4, eig = T)

p2 <- pcoa(dmr, design, 12, "Plot", "Stage", 2.4, kingdom) +
    theme(legend.position = "none")

############################### For fungal community
asv_file <- "../00.data/final_ASV/Fungi/ASV_raref.rds"
kingdom <- "Fungi"

asv <- readRDS(asv_file)
asv <- asv[, colnames(asv) %in% design$Sample_ID]

# Bray-Curtis distance
dis <- as.matrix(parDist(t(asv), method = "bray"))
dmr <- cmdscale(dis, k = 4, eig = T)

p3 <- pcoa(dmr, design, 12, "Plot", "Stage", 2.4, kingdom) +
    theme(legend.position = "none")

############################### For oomycetal community
asv_file <- "../00.data/final_ASV/Oomycetes/ASV_raref.rds"
kingdom <- "Oomycetes"

asv <- readRDS(asv_file)
asv <- asv[, colnames(asv) %in% design$Sample_ID]

# Bray-Curtis distance
dis <- as.matrix(parDist(t(asv), method = "bray"))
dmr <- cmdscale(dis, k = 4, eig = T)
## Compartment + Stage
p4 <- pcoa(dmr, design, 12, "Plot", "Stage", 2.4, kingdom) +
    theme(legend.position = "none")

## combine together
library(gridExtra)
all <- grid.arrange(p1, p_ph, p_clay,
                    p2, p3, p4,
                    nrow = 2, ncol = 3)
ggsave("Figure_S2.pdf", all, width = 12, height = 8)
