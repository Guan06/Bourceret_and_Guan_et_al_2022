#!/biodata/dep_psl/grp_rgo/tools/bin/Rscript
library(ggplot2)
library(ggpubr)
library(gridExtra)
library(agricolae)
source("../00.common_scripts/plot_settings.R")
source("./figure_s10_settings.R")
## for diversity
div_bac <- read.table("../00.data/alpha_diversity/alpha_Bacteria.txt",
                      header = T, sep = "\t")
div_fun <- read.table("../00.data/alpha_diversity/alpha_Fungi.txt",
                      header = T, sep = "\t")
div_oom <- read.table("../00.data/alpha_diversity/alpha_Oomycetes.txt",
                      header = T, sep = "\t")

design <- read.table("../00.data/design_48.txt", header = T, sep = "\t")

design <- design[, colnames(design) %in% c("Sample_ID", "Stage", "Management",
                                           "Genotype")]

div_bac <- div_bac[, colnames(div_bac) %in% c("Sample_ID", "Shannon")]
div_fun <- div_fun[, colnames(div_fun) %in% c("Sample_ID", "Shannon")]
div_oom <- div_oom[, colnames(div_oom) %in% c("Sample_ID", "Shannon")]

box_div <- function(x, design, prefix) {
    m <- merge(design, x)
    m$Man_Sta <- paste(m$Management, m$Stage, sep = "_")
    m$Man_Sta <- factor(m$Man_Sta, levels = order0, ordered = TRUE)

    m$Group <- paste(m$Stage, m$Management, m$Genotype, sep = "_")
    m$Group <- factor(m$Group, levels = order, ordered = TRUE)
    p <- box(m, "Genotype", "Man_Sta", "Group", "Shannon") +
        scale_x_discrete(limits = order) + ggtitle(prefix) +
        theme(plot.title = element_text(size = 12),
              legend.position = "none",
              axis.text.x = element_text(colour = "black", angle = 90,
              size = 8, hjust = 1))
    return(p)
}

p1 <- box_div(div_bac, design, "Bacteria")
p2 <- box_div(div_fun, design, "Fungi")
p3 <- box_div(div_oom, design, "Oomycetes")

all <- grid.arrange(p1, p2, p3, nrow = 1, ncol = 3)
ggsave("Figure_S10c.pdf", all, width = 10, height = 3)

shannon_sig <- function(x , design) {
    m <- merge(design, x)
    m$Man_Sta <- paste(m$Management, m$Stage, sep = "_")
    m$Group <- paste(m$Stage, m$Management, m$Genotype, sep = "_")

    sig1 <- box_sig(m, "Stage", "Shannon")
    sig2 <- box_sig(m, "Man_Sta", "Shannon")
    sig3 <- box_sig(m, "Group", "Shannon")

    sig <- rbind(sig1, sig2, sig3)
    sig$Sig <- ifelse(as.numeric(as.character(sig$Significance)) < 0.05,
                      TRUE, FALSE)
    return(sig)
}

sig_bac <- shannon_sig(div_bac, design)
sig_bac$Kingdom <- "Bacteria"
sig_fun <- shannon_sig(div_fun, design)
sig_fun$Kingdom <- "Fungi"
sig_oom <- shannon_sig(div_oom, design)
sig_oom$Kingdom <- "Oomycetes"

sig_tab <- rbind(sig_bac, sig_fun, sig_oom)
sig_tab$Significance <- as.numeric(as.character(sig_tab$Significance))
sig_tab$Sig <- ifelse(sig_tab$Significance < 0.05, TRUE, FALSE)
sig_tab$FDR <- p.adjust(sig_tab$Significance, method = "fdr")
sig_tab$Sig_FDR <- ifelse(as.numeric(as.character(sig_tab$FDR)) < 0.05,
                        TRUE, FALSE)

write.table(sig_tab, "Figure_S10c.txt", row.names = F, sep = "\t", quote = F)
