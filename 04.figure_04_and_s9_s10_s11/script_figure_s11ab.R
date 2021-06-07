#!/biodata/dep_psl/grp_rgo/guan/bin/Rscript

library(cowplot)

source("../00.common_scripts/plot_settings.R")

bac_file <- "../00.data/OTUs/bacteria/root_persistent_26.lst"
fun_file <- "../00.data/OTUs/fungi/root_persistent_27.lst"
bac_lst <- read.table(bac_file, header = F, sep = "\t")
fun_lst <- read.table(fun_file, header = F, sep = "\t")

###############################################################################
im_mse_bac <- read.table("Figure_S11_bac_mse.txt", header = T, sep = "\t")
im_mse_fun <- read.table("Figure_S11_fun_mse.txt", header = T, sep = "\t")

im_sum_bac <- data.frame(Taxon = rownames(im_mse_bac),
                         Total_MSE = rowSums(im_mse_bac))
im_sum_fun <- data.frame(Taxon = rownames(im_mse_fun),
                         Total_MSE = rowSums(im_mse_fun))

top_bac <- im_sum_bac[1:25, ]
top_bac$Group <- ifelse(top_bac$Taxon %in% bac_lst[, 1], "Stable", "Dynamic")

top_fun <- im_sum_fun[1:25, ]
top_fun$Group <- ifelse(top_fun$Taxon %in% fun_lst[, 1], "Stable", "Dynamic")

write.table(top_bac, "Figure_S11_bac_top.txt",
            quote = F, sep = "\t", row.names = F)
write.table(top_fun, "Figure_S11_fun_top.txt",
            quote = F, sep = "\t", row.names = F)

c_core <- c("Stable" = "salmon", "Dynamic" = "gray64")

### for plot a and b

p_a <- ggplot(top_bac, aes(Taxon, Total_MSE)) +
    geom_bar(aes(fill = Group), stat = 'identity') +
    scale_fill_manual(values = c_core) +
    scale_x_discrete(limits = rev(top_bac$Taxon)) +
    coord_flip() +
    main_theme +
    theme(legend.position = 'top') +
    labs(y = 'Total mean squred error')

p_b <- ggplot(top_fun, aes(Taxon, Total_MSE)) +
    geom_bar(aes(fill = Group), stat = 'identity') +
    scale_fill_manual(values = c_core) +
    scale_x_discrete(limits = rev(top_fun$Taxon)) +
    coord_flip() +
    main_theme +
    theme(legend.position = 'top') +
    labs(y = 'Total mean squared error')

p_ab <- plot_grid(p_a, p_b, nrow = 1)

ggsave("Figure_S11ab.pdf", p_ab, height = 5, width = 7)
