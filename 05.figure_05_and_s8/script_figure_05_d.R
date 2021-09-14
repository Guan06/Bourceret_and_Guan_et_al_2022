#!/biodata/dep_psl/grp_rgo/tools/bin/Rscript
library(ggplot2)
library(ggpubr)
library(gridExtra)
source("../00.common_scripts/plot_settings.R")
source("./figure_05_settings.R")

biomass <- read.table("../00.data/meta_data/biomass_noNA.txt",
                      header = T, sep = "\t")
design <- read.table("../00.data/design_48.txt", header = T, sep = "\t")

design$Man_Sta <- paste0(design$Management, "_", design$Stage)
design$Sta_Man_Gen <- paste0(design$Stage, "_",
                             design$Management, "_", design$Genotype)

design <- design[, colnames(design) %in% c("Sample_ID", "Stage", "Management",
                                           "Genotype",
                                           "Man_Sta", "Sta_Man_Gen")]

biomass <- biomass[, c(1, 2, 4, 7, 8)]
m <- merge(biomass, design)

lst_Stage <- as.character(unique(m$Stage))
lst_Management <- as.character(unique(m$Management))
lst_Genotype <- as.character(unique(m$Genotype))

for(n in 2 : 5){
    this_biomass <- m[, c(n, 6:10)]
    colnames(this_biomass)[6] <- "Group"
    this_biomass$Group <- factor(this_biomass$Group,
                                 levels = order,
                                 ordered = TRUE)
    this_biomass_name <- colnames(m)[n]
    p1 <- box(this_biomass, "Genotype", "Man_Sta", "Group", this_biomass_name) +
        scale_x_discrete(limits = order) +
        theme(legend.position = "none",
              axis.text.x = element_text(colour = "black", angle = 90,
                                         size = 8, hjust = 1))

    pdf <- paste0("Figure_5d_" , this_biomass_name, ".pdf")
    ggsave(pdf, p1, width = 3, height = 3.75)

    sig0 <- box_sig(this_biomass, "Genotype", this_biomass_name)
    sig1 <- c()
    sig2 <- c()
    sig3 <- c()
    for (s in lst_Stage) {
        this_stage <- this_biomass[this_biomass$Stage == s, ]
        this_sig1 <- box_sig(this_stage, "Management", this_biomass_name)
        this_sig1$Group1 <- paste0(s, "_", this_sig1$Group1)
        this_sig1$Group2 <- paste0(s, "_", this_sig1$Group2)
        sig1 <- rbind(sig1, this_sig1)
        for (mn in lst_Management) {
            this <- this_biomass[this_biomass$Stage == s, ]
            this <- this[this$Management == mn, ]
            this_sig <- box_sig(this, "Genotype", this_biomass_name)
            this_sig$Group1 <- paste0(s, "_", mn, "_", this_sig$Group1)
            this_sig$Group2 <- paste0(s, "_", mn, "_", this_sig$Group2)
            sig2 <- rbind(sig2, this_sig)
        }
        for (gt in lst_Genotype) {
            this <- this_biomass[this_biomass$Stage == s, ]
            this <- this[this$Genotype == gt, ]
            this_sig <- box_sig(this, "Management", this_biomass_name)
            this_sig$Group1 <- paste0(s, "_", gt, "_", this_sig$Group1)
            this_sig$Group2 <- paste0(s, "_", gt, "_", this_sig$Group2)
            sig3 <- rbind(sig3, this_sig)
        }
    }
    sig_all <- rbind(sig0, sig1, sig2, sig3)
    sig_all$Significance <- as.numeric(as.character(sig_all$Significance))
    sig_all$Sig <- ifelse(sig_all$Significance < 0.05, TRUE, FALSE)
    sig_all$FDR <- p.adjust(sig_all$Significance, method = "fdr")
    sig_all$Sig_FDR <- ifelse(as.numeric(as.character(sig_all$FDR)) < 0.05,
                    TRUE, FALSE)
    write.table(sig_all, paste0("Figure_5d_", this_biomass_name, ".txt"), quote = F,
                sep = "\t", row.names = F)
}
