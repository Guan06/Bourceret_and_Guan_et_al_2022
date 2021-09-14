#!/biodata/dep_psl/grp_rgo/guan/bin/Rscript

source("./function_get_network.R")
source("../00.common_scripts/plot_settings.R")

design_file <- "../00.data/design.txt"
design <- read.table(design_file, header = T, sep = "\t")
design <- design[design$Compartment == "Root", ]

bac_otu_file <- "../00.data/OTUs/bacteria/otu_table.rds"
fun_otu_file <- "../00.data/OTUs/fungi/otu_table.rds"
bac_otu <- readRDS(bac_otu_file)
fun_otu <- readRDS(fun_otu_file)

bac_otu <- apply(bac_otu, 2, function(x) x / sum(x))
fun_otu <- apply(fun_otu, 2, function(x) x / sum(x))

meta_file <- "../00.data/meta_data/metabolites_lipid.txt"
meta <- read.table(meta_file, header = T, sep="\t")
top_bac <- read.table("Figure_S7_bac_top_tax.txt", header = T, sep = "\t")
top_fun <- read.table("Figure_S7_fun_top_tax.txt", header = T, sep = "\t")
top_bac <- top_bac[, 1:6]
top_fun <- top_fun[, 1:6]

colnames(top_bac) <- colnames(top_fun) <- c("Taxon", "Kingdom", "Phylum",
                                            "Class", "Order", "Family")
###############################################################################
## compare the microbial RA and lipid change in VG and RP
lipid <- meta

design <- design[design$Genotype != "5_pht1;6", ]
samples <- intersect(lipid$Sample_ID, design$Sample_ID)
design <- design[design$Sample_ID %in% samples, ]

bac_otu <- apply(bac_otu, 2, function(x) x / sum(x))
fun_otu <- apply(fun_otu, 2, function(x) x / sum(x))

bac_otu2 <- bac_otu[rownames(bac_otu) %in% top_bac$Taxon,
                   colnames(bac_otu) %in% samples]
fun_otu2 <- fun_otu[rownames(fun_otu) %in% top_fun$Taxon,
                   colnames(fun_otu) %in% samples]

lipid <- lipid[lipid$Sample_ID %in% samples, ]
rownames(lipid) <- lipid$Sample_ID
lipid <- t(as.matrix(lipid[, -1]))

VG_samples <- design[design$Stage == "Vegetative", ]$Sample_ID
RP_samples <- design[design$Stage == "Reproductive", ]$Sample_ID

get_stage_mean <- function(x) {
    VG <- x[, colnames(x) %in% VG_samples]
    RP <- x[, colnames(x) %in% RP_samples]

    mean_df <- data.frame(Taxon = rownames(x),
                          VG_mean = rowSums(VG) / ncol(VG),
                          RP_mean = rowSums(RP) / ncol(RP))
}

mean_bac <- get_stage_mean(bac_otu2)
mean_fun <- get_stage_mean(fun_otu2)
mean_lipid <- get_stage_mean(lipid)

both <- rbind(mean_bac, mean_fun)
node_change <- rbind(both, mean_lipid)

node_change$Change <- ifelse(node_change$VG_mean > node_change$RP_mean,
                             "Decrease", "Increase")
###############################################################################

design <- read.table(design_file, header = T, sep = "\t")
design <- design[design$Compartment == "Root", ]

### for bacteria
bac_otu <- bac_otu[rownames(bac_otu) %in% top_bac$Taxon, ]
bac_otu <- bac_otu[, colSums(bac_otu) > 0]

bac_otu <- t(bac_otu)
bac_otu <- as.data.frame(bac_otu)
bac_otu$Sample_ID <- rownames(bac_otu)

bac_otu_meta <- merge(bac_otu, meta)

rownames(bac_otu_meta) <- bac_otu_meta$Sample_ID
bac_otu_meta <- bac_otu_meta[, -1]
bac_otu_meta <- as.matrix(bac_otu_meta[!is.na(rowSums(bac_otu_meta)), ])

this_sp <- rcorr(bac_otu_meta, type = "spearman")
net <- flattenCorrMatrix(this_sp$r, this_sp$P)
net[is.na(net)] <- 0
net <- net[net$edge != 0, ]
net <- net[net$p < 0.05, ]
net <- net[abs(net$edge) > 0.5, ]

node <- get_node_features(net)
node <- merge(node, top_bac)
node <- merge(node, node_change)

colnames(net)[1 : 3] <- c("source", "target", "weight")
colnames(node)[1] <- "id"
edge_file <- "Figure_S7c_bac_edge.txt"
node_file <- "Figure_S7c_bac_node.txt"
write.table(net, edge_file, quote = F, sep = "\t", row.names = F)
write.table(node, node_file, quote = F, sep = "\t", row.names = F)

### for fungi
fun_otu <- fun_otu[rownames(fun_otu) %in% top_fun$Taxon, ]
fun_otu <- fun_otu[, colSums(fun_otu) > 0]

fun_otu <- t(fun_otu)
fun_otu <- as.data.frame(fun_otu)
fun_otu$Sample_ID <- rownames(fun_otu)

fun_otu_meta <- merge(fun_otu, meta)

rownames(fun_otu_meta) <- fun_otu_meta$Sample_ID
fun_otu_meta <- fun_otu_meta[, -1]
fun_otu_meta <- as.matrix(fun_otu_meta[!is.na(rowSums(fun_otu_meta)), ])

this_sp <- rcorr(fun_otu_meta, type = "spearman")
net <- flattenCorrMatrix(this_sp$r, this_sp$P)
net[is.na(net)] <- 0
net <- net[net$edge != 0, ]
net <- net[net$p < 0.05, ]
net <- net[abs(net$edge) > 0.5, ]

node <- get_node_features(net)
node <- merge(node, top_fun)
node <- merge(node, node_change)

colnames(net)[1 : 3] <- c("source", "target", "weight")
colnames(node)[1] <- "id"
edge_file <- "Figure_S7c_fun_edge.txt"
node_file <- "Figure_S7c_fun_node.txt"
write.table(net, edge_file, quote = F, sep = "\t", row.names = F)
write.table(node, node_file, quote = F, sep = "\t", row.names = F)
