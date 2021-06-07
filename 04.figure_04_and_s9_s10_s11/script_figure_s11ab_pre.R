#!/biodata/dep_psl/grp_rgo/guan/bin/Rscript

source("../00.common_scripts/plot_settings.R")
source("./function_get_regressor.R")
design_file <- "../00.data/design.txt"

bac_otu_file <- "../00.data/OTUs/bacteria/otu_table.rds"
fun_otu_file <- "../00.data/OTUs/fungi/otu_table.rds"
bac_file <- "../00.data/OTUs/bacteria/root_persistent_26.lst"
fun_file <- "../00.data/OTUs/fungi/root_persistent_27.lst"
meta_file <- "../00.data/meta_data/metabolites_lipid.txt"

design <- read.table(design_file, header = T, sep = "\t")
bac_lst <- read.table(bac_file, header = F, sep = "\t")
bac_otu <- readRDS(bac_otu_file)
fun_lst <- read.table(fun_file, header = F, sep = "\t")
fun_otu <- readRDS(fun_otu_file)

bac_otu <- apply(bac_otu, 2, function(x) x / sum(x))
fun_otu <- apply(fun_otu, 2, function(x) x / sum(x))

design <- design[design$Compartment == "Root", ]
meta <- read.table(meta_file, header = T, sep="\t")

## regressor by Random Forest, will take some time
im_mse_bac <- get_regressor(bac_otu, meta, "Figure_S11_bac")
im_mse_bac <- im_mse_bac[order(rowSums(im_mse_bac),
                               decreasing = TRUE), ]

write.table(im_mse_bac, "Figure_S11_bac_mse.txt",
            quote = F, sep = "\t")

## the regressor will take long time to run
im_mse_fun <- get_regressor(fun_otu, meta, "Figure_S11_fun")
im_mse_fun <- im_mse_fun[order(rowSums(im_mse_fun),
                               decreasing = TRUE), ]

write.table(im_mse_fun, "Figure_S11_fun_mse.txt",
            quote = F, sep = "\t")
