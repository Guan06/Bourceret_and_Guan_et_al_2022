#!/biodata/dep_psl/grp_rgo/tools/bin/Rscript
library(ggplot2)
library(gridExtra)
source("../00.common_scripts/plot_settings.R")

design <- read.table("../00.data/design_48.txt",
                     header = T, sep = "\t")
design <- design[, colnames(design) %in% c("Sample_ID", "Stage", "Management",
                                           "Genotype")]

biomass <- read.table("../00.data/meta_data/biomass_noNA.txt",
                      header = T, sep = "\t")
sugar <- read.table("../00.data/meta_data/metabolites_sugar.txt",
                    header = T, sep = "\t")
aa <- read.table("../00.data/meta_data/metabolites_aa.txt",
                 header = T, sep = "\t")
lipid <- read.table("../00.data/meta_data/metabolites_lipid.txt",
                    header = T, sep = "\t")
ionomic <- read.table("../00.data/meta_data/ionomic_20.txt",
                      header = T, sep = "\t")

lst_Stage <- c("VG", "RP")
lst_Man <- c("NK", "NPK")

## filter biomass
biomass <- biomass[, colnames(biomass) %in% c("Sample_ID", "dry.weight.cob",
                                              "dry.weight.shoot", "stem.length",
                                              "total.leaf.number")]

get_sig <- function(x) {
    m <- merge(design, x)
    nc0 <- ncol(design) + 1
    nc <- ncol(m)
    out <- c()

    for (n in nc0 : nc) {
        this <- m[, c(1 : 5, n)]
        this_name <- colnames(m)[n]

        for (s in lst_Stage) {
            this_stage <- this[this$Stage == s, ]
            for (mn in lst_Man) {
                this_mn <- this_stage[this_stage$Management == mn, ]
                if (sum(this_mn[[this_name]]) == 0) next
                ## get significance and fold change of mean
                wild <- this_mn[this_mn$Genotype == "B73", ]
                mutant <- this_mn[this_mn$Genotype == "pht1;6", ]

                p_value <- wilcox.test(wild[[this_name]], mutant[[this_name]],
                                       alternative = "two.sided")$p.value
                if (p_value > 0.05) next
                if (is.na(p_value)) next
                fc <- mean(wild[[this_name]]) / mean(mutant[[this_name]])

                g1 <- paste0(s, "_", mn, "_", "B73")
                g2 <- paste0(s, "_", mn, "_", "pht1;6")

                this_out <- c(g1, g2, p_value, fc, this_name)
                out <- rbind(out, this_out)
            }
        }
    }
    colnames(out) <- c("G1", "G2", "p_value", "FC", "Meta")
    return(out)
}

sig_biomass <- get_sig(biomass)
sig_sugar <- get_sig(sugar)
sig_aa <- get_sig(aa)
sig_lipid <- get_sig(lipid)
sig_ionomic <- get_sig(ionomic)

lst_biomass <- unique(sig_biomass[, 5])
lst_sugar <- unique(sig_sugar[, 5])
lst_aa <- unique(sig_aa[, 5])
lst_ionomic <- unique(sig_ionomic[, 5])
lst_lipid <- unique(sig_lipid[, 5])

## get subset of meta data which is at least significant for one comparison
## and calculate the corresponding FC
get_sub <- function(x, lst) {
    m <- merge(design, x)
    sub <- cbind(m[, 1:5], m[, colnames(m) %in% lst])
    nc0 <- ncol(design) + 1
    nc <- ncol(sub)
    out <- c()

    for (s in lst_Stage) {
        this_stage <- sub[sub$Stage == s, ]
        for (mn in lst_Man) {
            this_mn <- this_stage[this_stage$Management == mn, ]
            for (n in nc0 : nc) {
                this_meta <- this_mn[, c(1 : 5, n)]
                this_name <- colnames(this_mn)[n]

                if (sum(this_meta[[this_name]]) == 0) {
                    fc <- 0
                    p_value <- 1
                } else {
                    wild <- this_meta[this_meta$Genotype == "B73", ]
                    wild_mean <- mean(wild[[this_name]])

                    mutant <- this_meta[this_meta$Genotype == "pht1;6", ]
                    mutant_mean <- mean(mutant[[this_name]])

                    if (mutant_mean < wild_mean) {
                        fc <- - (wild_mean / mutant_mean)
                    } else {
                        fc <- mutant_mean / wild_mean
                    }
                    p_value <- wilcox.test(wild[[this_name]], mutant[[this_name]],
                                           alternative = "two.sided")$p.value
                }
                this_out <- c(s, mn, "B73", "pht1;6", p_value, fc, this_name)
                out <- rbind(out, this_out)
            }
        }
    }
    colnames(out) <- c("Stage", "Management", "Genotype_1", "Genotype_2",
                       "p_value", "Fold_change", "Meta")
    return(out)
}

sub_biomass <- get_sub(biomass, lst_biomass)
sub_sugar <- get_sub(sugar, lst_sugar)
sub_aa <- get_sub(aa, lst_aa)
sub_ionomic <- get_sub(ionomic, lst_ionomic)
sub_lipid <- get_sub(lipid, lst_lipid)

sub_all <- as.data.frame(rbind(sub_biomass, sub_sugar, sub_aa, sub_ionomic))
sub_tab <- rbind(sub_all, sub_lipid)
write.table(sub_tab, "Figure_4a.txt", quote = F, sep = "\t", row.names = F)

nk <- sub_all[sub_all$Management == "NK", ]
npk <- sub_all[sub_all$Management == "NPK", ]

nk$Fold_change <- as.numeric(as.character(nk$Fold_change))

p_nk <- ggplot(nk, aes(x = Meta, y = Fold_change, fill = Stage)) +
    geom_bar(width = 0.6, color = c_grey,
             position = position_dodge(width = 0.6),
             stat = "identity") +
    scale_x_discrete(limits = rev(unique(nk$Meta))) +
    coord_flip() + scale_fill_manual(values =c(c_grey, "white")) +
    main_theme + theme(legend.position = "top",
                       axis.text.y = element_text(size = 8),
                       axis.line.y = element_line(color = NA)) +
    geom_hline(yintercept = 0) +
    geom_hline(yintercept = -1, linetype = "dashed", color = "gold") +
    geom_hline(yintercept = 1, linetype = "dashed", color = "gold")


npk$Fold_change <- as.numeric(as.character(npk$Fold_change))

p_npk <- ggplot(npk, aes(x = Meta, y = Fold_change, fill = Stage)) +
    geom_bar(width = 0.6, color = c_black,
             position = position_dodge(width = 0.6), stat = "identity") +
    scale_x_discrete(limits = rev(unique(nk$Meta))) +
    coord_flip() + scale_fill_manual(values =c(c_black, "white")) +
    main_theme + theme(legend.position = "top",
                       axis.text.y = element_text(size = 8),
                       axis.line.y = element_line(color = NA)) +
    geom_hline(yintercept = 0) +
    geom_hline(yintercept = -1, linetype = "dashed", color = "gold") +
    geom_hline(yintercept = 1, linetype = "dashed", color = "gold")

## combine together
library(gridExtra)
all <- grid.arrange(p_nk, p_npk, nrow = 1, ncol = 2)
ggsave("Figure_4a.pdf", all, width = 10, height = 7)
