#!/biodata/dep_psl/grp_rgo/tools/bin/Rscript

library(ggplot2)

source("../00.common_scripts/plot_settings.R")

adonis <- read.table("../00.data/adonis.txt",
                     header = T, sep = "\t")

###############################################################################
## For Bacteria
adonis$Factor <- factor(adonis$Factor, levels = adonis$Factor)

p1 <- ggplot(adonis, aes(x = Factor, y = Bacteria)) + main_theme +
    geom_bar(stat = "identity", fill = "gray61") +
    theme(axis.text.x = element_text(angle = 90, size = 6,
                                     vjust = 0.5, hjust = 1)) +
    labs(x = "", y = "") + ggtitle("Bacteria")

## For Fungi
fp1 <- ggplot(adonis, aes(x = Factor, y = Fungi)) + main_theme +
    geom_bar(stat = "identity", fill = "gray61") +
    theme(axis.text.x = element_text(angle = 90, size = 6,
                                     vjust = 0.5, hjust = 1)) +
    labs(x = "", y = "") +ggtitle("Fungi")

## For Oomycetes
op1 <- ggplot(adonis, aes(x = Factor, y = Oomycetes)) + main_theme +
    geom_bar(stat = "identity", fill = "gray61") +
    theme(axis.text.x = element_text(angle = 90, size = 6,
                                     vjust = 0.5, hjust = 1)) +
    labs(x = "", y = "") + ggtitle("Oomycetes")

## output
library(gridExtra)

all_p1 <- grid.arrange(p1, fp1, op1, nrow = 1, ncol = 3)
ggsave("Figure_2b.pdf", all_p1, width = 10, height = 3)
