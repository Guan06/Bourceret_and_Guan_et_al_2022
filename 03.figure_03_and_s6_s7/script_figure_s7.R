#!/biodata/dep_psl/grp_rgo/tools/bin/Rscript
library(ggplot2)
library(ggpubr)
library(gdata)

source("../00.common_scripts/plot_settings.R")
source("../00.common_scripts/box_plotting.R")

myc_file <- "../00.data/meta_data/colonization_myc.txt"
design_file <- "../00.data/design.txt"
asv_file <- "../00.data/final_ASV/Fungi/ASV_table_1104.rds"
tax_file <- "../00.data/final_ASV/Fungi/ASV_taxonomy.txt"

asv <- readRDS(asv_file)
asv_norm <- apply(asv, 2, function(x) x /sum(x))
tax <- read.table(tax_file, header = T, sep = "\t")

## for myc
order_myc <- c("Vegetative_NK_1_B73",
                  "Vegetative_NPK_1_B73",
                  "Vegetative_CONMIN_1_B73",
                  "Vegetative_BIODYN_1_B73",
                  "Vegetative_NK_5_pht1;6",
                  "Vegetative_NPK_5_pht1;6",
                  "Vegetative_CONMIN_5_pht1;6",
                  "Vegetative_BIODYN_5_pht1;6",
                  "Reproductive_NK_1_B73",
                  "Reproductive_NPK_1_B73",
                  "Reproductive_CONMIN_1_B73",
                  "Reproductive_BIODYN_1_B73")

myc <- read.table(myc_file, header = T, sep = "\t")
design <- read.table(design_file, header = T, sep = "\t")

design <- design[design$Sample_ID %in% myc$Sample_ID, ]
design$Stage_Management_Genotype <- paste(design$Stage, design$Management,
                                          design$Genotype, sep = "_")
design <- design[, colnames(design) %in% c("Sample_ID", "Stage", "Management",
                                           "Genotype", "Stage_Management_Genotype")]
## for p__Glomeromycota relative abundance
tax$ASV_ID <- rownames(tax)
tax <- tax[tax$Phylum == "p__Glomeromycota", ]

asv_norm_glo <- asv_norm[rownames(asv_norm) %in% tax$ASV_ID, ]
glo_ra <- data.frame(Sample_ID = colnames(asv_norm_glo),
                     RA = colSums(asv_norm_glo))

ra <- merge(glo_ra, design)

library(dplyr)
ra_sum <- ra %>% group_by(Stage_Management_Genotype) %>%
    summarise(mean = mean(RA), sd = sd(RA))

a <- ggplot(ra_sum) + geom_bar(aes(x = Stage_Management_Genotype, y = mean),
                              stat = "identity") + main_theme +
        scale_x_discrete(limits = order_myc) +
        theme(legend.position = "none",
              axis.text.x = element_blank())

m <- merge(design, myc)
m <- drop.levels(m)
m <- m[, c(1:5, 7)]

b <- box(m, "Management", "Stage", "Stage_Management_Genotype",
              "myc_colonization_degree") +
        scale_x_discrete(limits = order_myc) +
        scale_shape_manual(values = c(1, 16)) +
        theme(legend.position = "none",
        axis.text.x = element_blank())

## combine
library(gridExtra)

s7 <- grid.arrange(a, b, nrow = 2, ncol = 1)
ggsave("Figure_S7.pdf", s7)

lst_Stage <- as.character(unique(m$Stage))
lst_Management <- as.character(unique(m$Management))
lst_Genotype <- as.character(unique(m$Genotype))

## Test for a ################################################################
sig0 <- box_sig(ra, "Genotype", "RA")

sig1 <- c()
s <- "Vegetative"
this_stage <- ra[ra$Stage == s, ]
sig1 <- box_sig(this_stage, "Genotype", "RA")
sig1$Group1 <- paste0(s, "_", sig1$Group1)
sig1$Group2 <- paste0(s, "_", sig1$Group2)

sig2 <- c()
for (s in lst_Stage) {
        ## only B73 for reproductive stage
    this_stage <- ra[ra$Stage == s, ]
    ## for management test
    for (g in lst_Genotype) {
        this_genotype <- this_stage[this_stage$Genotype == g, ]
        if (nrow(this_genotype) == 0) next
        this_sig2 <- box_sig(this_genotype, "Management", "RA")
        this_sig2$Group1 <- paste0(s, "_", g, "_", this_sig2$Group1)
        this_sig2$Group2 <- paste0(s, "_", g, "_", this_sig2$Group2)
        sig2 <- rbind(sig2, this_sig2)
    }
}
## test B73 for VG and RP
sig3 <- c()
b73 <- ra[ra$Genotype == "1_B73", ]
sig3 <- box_sig(b73, "Stage", "RA")
sig3$Group1 <- paste0("1_B73_", sig3$Group1)
sig3$Group2 <- paste0("1_B73_", sig3$Group2)

sig4 <- c()
for (mn in lst_Management) {
    this_man <- b73[b73$Management == mn, ]
    this_sig4 <- box_sig(this_man, "Stage", "RA")
    this_sig4$Group1 <- paste0("1_B73_", mn, "_", this_sig4$Group1)
    this_sig4$Group2 <- paste0("1_B73_", mn, "_", this_sig4$Group2)
    sig4 <- rbind(sig4, this_sig4)
}

sig_a <- rbind(sig0, sig1, sig3, sig4, sig2)
sig_a$Significance <- as.numeric(as.character(sig_a$Significance))
sig_a$Sig <- ifelse(sig_a$Significance < 0.05, TRUE, FALSE)
sig_a$FDR <- p.adjust(sig_a$Significance, method = "fdr")
sig_a$Sig_FDR <- ifelse(as.numeric(as.character(sig_a$FDR)) < 0.05,
                        TRUE, FALSE)

write.table(sig_a, "Figure_S7a.txt", quote = F, sep = "\t", row.names = F)

## Test for b ################################################################
sig0 <- box_sig(m, "Genotype", "myc_colonization_degree")

sig1 <- c()
s <- "Vegetative"
this_stage <- m[m$Stage == s, ]
sig1 <- box_sig(this_stage, "Genotype", "myc_colonization_degree")
sig1$Group1 <- paste0(s, "_", sig1$Group1)
sig1$Group2 <- paste0(s, "_", sig1$Group2)

sig2 <- c()
for (s in lst_Stage) {
    ## only B73 for reproductive stage
    this_stage <- m[m$Stage == s, ]
    ## for management test
    for (g in lst_Genotype) {
        this_genotype <- this_stage[this_stage$Genotype == g, ]
        if (nrow(this_genotype) == 0) next
        this_sig2 <- box_sig(this_genotype, "Management", "myc_colonization_degree")
        this_sig2$Group1 <- paste0(s, "_", g, "_", this_sig2$Group1)
        this_sig2$Group2 <- paste0(s, "_", g, "_", this_sig2$Group2)
        sig2 <- rbind(sig2, this_sig2)
    }
}

## test B73 for VG and RP
sig3 <- c()
b73 <- m[m$Genotype == "1_B73", ]
sig3 <- box_sig(b73, "Stage", "myc_colonization_degree")
sig3$Group1 <- paste0("1_B73_", sig3$Group1)
sig3$Group2 <- paste0("1_B73_", sig3$Group2)

# for each management
sig4 <- c()
for (mn in lst_Management) {
    this_man <- b73[b73$Management == mn, ]
    this_sig4 <- box_sig(this_man, "Stage", "myc_colonization_degree")
    this_sig4$Group1 <- paste0("1_B73_", mn, "_", this_sig4$Group1)
    this_sig4$Group2 <- paste0("1_B73_", mn, "_", this_sig4$Group2)
    sig4 <- rbind(sig4, this_sig4)
}

sig_b <- rbind(sig0, sig1, sig3, sig4, sig2)
sig_b$Significance <- as.numeric(as.character(sig_b$Significance))
sig_b$Sig <- ifelse(sig_b$Significance < 0.05, TRUE, FALSE)
sig_b$FDR <- p.adjust(sig_b$Significance, method = "fdr")
sig_b$Sig_FDR <- ifelse(as.numeric(as.character(sig_b$FDR)) < 0.05,
                        TRUE, FALSE)
write.table(sig_b,"Figure_S7b.txt", quote = F, sep = "\t", row.names = F)
