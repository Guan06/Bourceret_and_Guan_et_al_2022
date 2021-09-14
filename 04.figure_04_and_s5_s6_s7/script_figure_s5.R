#!/biodata/dep_psl/grp_rgo/guan/bin/Rscript

library(pheatmap)
library(dplyr)
source("../00.common_scripts/plot_settings.R")

design_file <- "../00.data/design.txt"
response_file <- "../00.data/meta_data/metabolites_lipid.txt"

design <- read.table(design_file, header = T, sep = "\t")
design <- design[design$Genotype != "5_pht1;6", ]

des <- design[, c("Sample_ID", "Stage", "Management", "Genotype")]
des$Group <- paste0(des$Stage, "_", des$Management, "_", des$Genotype)
rownames(des) <- des$Sample_ID

response <- read.table(response_file, header = T, sep="\t")
rownames(response) <- response$Sample_ID
response <- response[, -1]

## for lipid, remove sample with NA
response <- response[!(is.na(rowSums(response))), ]

des <- des[des$Sample_ID %in% rownames(response), ]
response <- response[match(rownames(des), rownames(response)), ]
response$Group <- des$Group

res_mean <- as.data.frame(response %>% group_by(Group) %>%
    summarise(across(everything(), list(mean))))

rownames(res_mean) <- res_mean$Group
res_mean <- res_mean[, -1]
res_mean <- log(res_mean)

## range normalize the ionomic data
get_range <- function(x) {
    (x - min(x)) / (max(x) - min(x))
}

res2 <- apply(res_mean, 2, get_range)
res2 <- t(res2)

## for each genotype
des2 <- unique(des[, c("Group", "Genotype", "Management", "Stage")])
rownames(des2) <- des2$Group
des2 <- des2[, -1]
greens <- brewer.pal(n = 9, name = "Greens")
c_Gen <- c("1_B73" = set2[3], "2_DK105" = set3[3], "3_PH207" = set2[5],
           "4_F2" = set3[4])
my_color = list(
    Genotype = c_Gen,
    Management = c_Man,
    Stage = c("Vegetative" = greens[3], "Reproductive" = greens[6])
)

rownames(res2) <- unlist(strsplit(rownames(res2), "_"))[seq(1,by = 2,
                                                            len = nrow(res2))]
## Sanzo Wasa 004
## https://sanzo-wada.dmbk.io/combination/004 
c_isabella <- "#c3a55c"
c_red_violet <- "#3c00a3"
cp <- colorRampPalette(c(c_red_violet, "white", c_isabella))(50)
p_all <- pheatmap(res2,
#                  color = cp,
                  cluster_rows = FALSE, show_colnames = FALSE,
                  annotation_colors = my_color,
                  annotation_col = des2, fontsize = 5.6)

ggsave("Figure_S5.pdf", p_all, width = 6, height = 5)
