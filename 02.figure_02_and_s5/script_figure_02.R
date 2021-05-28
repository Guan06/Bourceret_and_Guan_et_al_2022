#!/biodata/dep_psl/grp_rgo/tools/bin/Rscript

library(ggplot2)
library(cowplot)
source("../00.common_scripts/plot_settings.R")
source("../00.common_scripts/violin_plotting.R")

bac_lst <- "../00.data/OTUs/bacteria/widespread_persistent_15.lst"
bac_otu <- "../00.data/OTUs/bacteria/otu_table.rds"

fun_lst <- "../00.data/OTUs/fungi/widespread_persistent_24.lst"
fun_otu <- "../00.data/OTUs/fungi/otu_table.rds"

design_file <- "../00.data/design.txt"
design <- read.table(design_file, header = T, sep = "\t")

a_Sta <- c("Vegetative" = 0.25, "Reproductive" = 0.5)

get_ra_df <- function(x, lst, des = design) {
    otu <- apply(x, 2, function(x) x / sum(x))
    otu <- otu[rownames(otu) %in% lst, ]
    RA_DF <- c()
    for (com in c("Rhizosphere", "Root")) {
        for (st in c("Vegetative", "Reproductive")) {
            this_des <- des[(des$Compartment == com) &
                            (des$Stage == st), ]
            this_otu <- otu[, colnames(otu) %in% this_des$Sample_ID]
            print(paste(com, st, ncol(this_otu), " samples."))
            ra <- sum(colSums(this_otu)) / ncol(this_otu)
            print(ra)

            ra_df <- data.frame(Sample_ID = colnames(this_otu),
                            RA = colSums(this_otu))

            ra_df <- merge(ra_df, this_des)

            ra_df <- ra_df[, c("Sample_ID", "RA", "Compartment", "Stage", "Management",
                           "Compartment_Stage")]

            RA_DF <- rbind(RA_DF, ra_df)
        }
    }

    p <- ggplot(RA_DF, aes(Compartment_Stage, y = RA)) +
        geom_violin(aes(fill = Compartment, alpha = Stage), color = "NA") +
        geom_jitter(aes(color = Management, shape = Stage),
                    size = 1.8) +
        scale_fill_manual(values = c_Com) +
        scale_alpha_manual(values = a_Sta) +
        scale_colour_manual(values = c_Man) +
        scale_shape_manual(values = s_Sta) +
        main_theme +
        labs(x = "", y = "Aggregated relative abundance") +
        theme(axis.text.y = element_text(size = 18),
              axis.title.y = element_text(size = 18),
              axis.text.x = element_text(angle = 90, vjust = 0.5,
                                         hjust = 1, size = 18),
              legend.position = "none") +
        scale_x_discrete(limits = c("Rhizosphere_Vegetative",
                                    "Rhizosphere_Reproductive",
                                    "Root_Vegetative",
                                    "Root_Reproductive"),
                         labels = c("Rhizos_VG", "Rhizos_RP", "Root_VG", "Root_RP"))
    return(p)
}

bac <- readRDS(bac_otu)
fun <- readRDS(fun_otu)
lst1 <- read.table(bac_lst)
lst2 <- read.table(fun_lst)

lst1 <- lst1$V1
lst2 <- lst2$V1

p_bac <- get_ra_df(bac, lst1)
p_fun <- get_ra_df(fun, lst2)

p <- plot_grid(p_bac, p_fun, nrow = 2, ncol = 1, align = "v", axis = "l")
ggsave("Figure_2.pdf", p, height = 12, width = 6)
