#!/biodata/dep_psl/grp_rgo/tools/bin/Rscript

library(plyr)
library(parallelDist)
library(dplyr)
library(tidyr)
library(gridExtra)

source("../00.common_scripts/plot_settings.R")
source("../00.common_scripts/pcoa_plotting.R")
source("../00.common_scripts/hist_settings.R")

asv_file <- "../00.data/final_ASV/Fungi/ASV_table_1104.rds"
design_file <- "../00.data/design.txt"
tax_file <- "../00.data/final_ASV/Fungi/ASV_taxonomy.txt"
alpha_file <- "../00.data/alpha_diversity/alpha_Fungi.txt"
kingdom <- "Fungi"

###############################################################################
## For relative abundance histogram
asv <- readRDS(asv_file)
design <- read.table(design_file, header = T, sep = "\t")
tax <- read.table(tax_file, header = T, sep = "\t")

# Get the intersection of asv and design file
inter <- intersect(as.character(colnames(asv)), as.character(design$Sample_ID))
asv <- asv[, colnames(asv) %in% inter]
design <- design[design$Sample_ID %in% inter, ]

asv <- asv[, match(design$Sample_ID, colnames(asv))]

# normalize by colSums
asv <- apply(asv, 2, function(x) x / sum(x))

asv_tab <- as.table(as.matrix(asv))
asv_tab <- prop.table(asv_tab, 2)

asv_tab <- asv_tab * 1000 / 6
asv_tab <- as.data.frame(asv_tab)
colnames(asv_tab) <- c("ASV_ID", "Sample_ID", "Freq")

tax$ASV_ID <- rownames(tax)

asv_tax <- join(asv_tab, tax)
# dim = 10907520 * 10

asv_tax$Phylum <- as.character(asv_tax$Phylum)

## Use Unassigned for unassigned ASVs
idx3 <- is.na(asv_tax$Phylum)
asv_tax$Phylum[idx3] <- "Unassigned"

## Use others for phyla not included
idx4 <- !(asv_tax$Phylum %in% Fun_phyla)
asv_tax$Phylum[idx4] <- "Others"

asv_tax <- asv_tax[, c(1 : 3, 5)]
asv_tax <- asv_tax[asv_tax$Phylum %in% Fun_phyla_2, ]
# dim = 10907520 * 4

asv_tax_sum <- aggregate(asv_tax$Freq,
                         by = list(Phylum = asv_tax$Phylum,
                                   Sample_ID = asv_tax$Sample_ID),
                         FUN = sum)

## join with design file
col2 <- which(colnames(design) == "Compartment_Stage_Management_Plot_Genotype")
design <- design[, c(1, col2)]
colnames(design)[2] <- "Group2"

asv_tax_sum <- join(asv_tax_sum, design)

asv_tax_sum$Group2 <- factor(asv_tax_sum$Group2,
                             levels = order_new,
                             ordered = TRUE)

### Merge duplicates together

p <- ggplot(data = asv_tax_sum,
             aes(x = Group2, y = x,
                 fill = Phylum)) +
    geom_bar(stat = "identity", position = "stack") +
    scale_x_discrete(limits = order_new) +
    main_theme + scale_fill_manual(values = as.character(c_Fun$color)) +
    theme(legend.position = "top",
          axis.text.x = element_text(colour = "black", angle = 90,
                                     size = 4, hjust = 1)) +
    labs(x = "", y = "")

ggsave("Figure_S6a.pdf", p, width = 16, height = 10)

###############################################################################
## for B and C, phylum based PCoA and density plot
asv <- readRDS(asv_file)
design <- read.table(design_file, header = T, sep = "\t")
tax <- read.table(tax_file, header = T, sep = "\t")

# Tidy up tax to make it the same as histogram
tax$Phylum <- as.character(tax$Phylum)
tax$Phylum[is.na(tax$Phylum)] <- "Unassigned"

tax <- tax[, 1 : 2]
tax <- tax[rownames(tax) %in% rownames(asv), ]
# dim = 9880 * 2

group <- as.factor(tax$Phylum)
group_lst <- as.character(levels(group))
asv <- asv[match(rownames(tax), rownames(asv)), ]
asv <- apply(asv, 2, function(x) x / sum(x))
asv_tax <- apply(asv, 2, function(x) rowsum(x, group))
rownames(asv_tax) <- group_lst

bc <- parDist(t(asv_tax), method = "bray")
bc_mat <- as.matrix(bc)
dmr <- cmdscale(bc_mat, k = 4, eig = T)

design <- design[match(rownames(bc_mat), design$Sample_ID), ]
p_s6b <- pcoa(dmr, design, 12, "Compartment", "Stage", 1.2, kingdom) +
    theme(legend.position = "none")

## for density plot
dis_asv_mat <- as.matrix(parDist(t(asv), method = "bray"))
dis_tax_mat <- bc_mat

## filter design
design <- design[, colnames(design) %in% c("Sample_ID", "Compartment", "Stage",
                                            "Compartment_Stage")]

design_asv <- design[, colnames(design) %in% c("Sample_ID",
                                               "Compartment_Stage")]
colnames(design_asv) <- c("Compare", "Compare_Group")

############################################  at the ASV level
asv_dis <- c()
for ( i in 1: nrow(dis_asv_mat)){
    this <- data.frame(Sample_ID = rownames(dis_asv_mat)[i],
                       Compare = colnames(dis_asv_mat),
                       Dis = dis_asv_mat[i, ])
    this <- merge(this, design_asv)
    this <- this %>% group_by(Compare_Group) %>%
            mutate(Dis_Mean = mean(Dis))

    this_mean <- this[, colnames(this) %in%
            c("Sample_ID", "Compare_Group", "Dis_Mean")]
    this_mean <- unique(this_mean)

    asv_dis <- rbind(asv_dis, this_mean)
}

this_asv_filter <- asv_dis[asv_dis$Compare_Group == "Bulksoil_before_sowing" |
                           asv_dis$Compare_Group == "Root_Reproductive", ]

this_asv_filter <- this_asv_filter %>% spread(Compare_Group, Dis_Mean)
this_asv_filter <- merge(this_asv_filter, design)
this_asv_filter$Compartment <- factor(this_asv_filter$Compartment,
                                      levels = c_Com$group)
this_asv_filter$Stage <- factor(this_asv_filter$Stage,
                                levels = s_Sta$group)
this_asv_filter$Compartment_Stage <- factor(this_asv_filter$Compartment_Stage)

############################################  at the Phylum level
tax_dis <- c()
for ( i in 1: nrow(dis_tax_mat)){
    this <- data.frame(Sample_ID = rownames(dis_tax_mat)[i],
                       Compare = colnames(dis_tax_mat),
                       Dis = dis_tax_mat[i, ])
    this <- merge(this, design_asv)
    this <- this %>% group_by(Compare_Group) %>%
    mutate(Dis_Mean = mean(Dis))

    this_mean <- this[, colnames(this) %in%
            c("Sample_ID", "Compare_Group", "Dis_Mean")]
    this_mean <- unique(this_mean)

    tax_dis <- rbind(tax_dis, this_mean)
}

this_tax_filter <- tax_dis[tax_dis$Compare_Group == "Bulksoil_before_sowing" |
                           tax_dis$Compare_Group == "Root_Reproductive", ]
this_tax_filter <- this_tax_filter %>% spread(Compare_Group, Dis_Mean)
this_tax_filter <- merge(this_tax_filter, design)
this_tax_filter$Compartment <- factor(this_tax_filter$Compartment,
                                      levels = c_Com$group)
this_tax_filter$Stage <- factor(this_tax_filter$Stage,
                                levels = s_Sta$group)
this_tax_filter$Compartment_Stage <- factor(this_tax_filter$Compartment_Stage)

##########################################  density plot
greens <- brewer.pal(n = 9, name = "Greens")
oranges <- brewer.pal(n = 9, name = "Oranges")
br <- brewer.pal(n = 11, name = "BrBG")

############## at the ASV level
## reduce 2D to 1D by calculating Euclidean distance
nr <- nrow(this_asv_filter)
ASV <- c()
## set sample 558H (Bulksoil_before_sowing sample) as control
min_x <- min(this_asv_filter[, 2])
x0 <- min_x
y0 <- this_asv_filter[this_asv_filter[, 2] == min_x, 3]

for (n in 1 : nr ){
    ## for Bulksoil_before_sowing
    xn <- this_asv_filter[n, 2]
    ## for Root_Reproductive
    yn <- this_asv_filter[n, 3]
    dis <- sqrt((xn - x0) ** 2 + (yn - y0) **2)
    this <- c(as.character(this_asv_filter[n, 1]),
              as.character(this_asv_filter[n, 6]), dis)
    ASV <- rbind(ASV, this)
}

colnames(ASV) <- c("Sample_ID", "Compartment_Stage", "Distance")
ASV <- as.data.frame(ASV)
ASV$Distance <- as.numeric(as.character(ASV$Distance))
ASV$Compartment_Stage <- as.character(ASV$Compartment_Stage)

df_mean <- ASV %>% group_by(Compartment_Stage) %>%
    summarise(Mean = mean(Distance))
write.table(df_mean, "Figure_S6c_ASV_mean.txt",
            quote = F, sep = "\t", row.names = F)

lines_mean <- as.numeric(as.character(df_mean$Mean))
colors <- c(br[c(1, 2, 4)], oranges[c(6, 3)], greens[c(6, 3)])

p_s6c_1 <- ggplot(ASV, aes(Distance, fill = Compartment_Stage)) +
        geom_density(alpha = 0.7, size = 0.2, color = "wheat4") +
        scale_fill_manual(values = colors) +
        main_theme +
        theme(legend.position = "none", axis.title = element_blank()) +
        geom_vline(xintercept = lines_mean, color = colors,
                   linetype = "longdash")

############## at the phylum level
nr <- nrow(this_tax_filter)
TAX <- c()

min_x <- min(this_tax_filter[, 2])
x0 <- min_x
y0 <- this_asv_filter[this_tax_filter[, 2] == min_x, 3]

for( n in 1 : nr  ){
    xn <- this_tax_filter[n, 2]
    yn <- this_tax_filter[n, 3]
    dis <- sqrt((xn - x0) ** 2 + (yn - y0) **2)
    this <- c(as.character(this_tax_filter[n, 1]),
              as.character(this_tax_filter[n, 6]), dis)
    TAX <- rbind(TAX, this)
}

colnames(TAX) <- c("Sample_ID", "Compartment_Stage", "Distance")
TAX <- as.data.frame(TAX)
TAX$Distance <- as.numeric(as.character(TAX$Distance))
TAX$Compartment_Stage <- as.character(TAX$Compartment_Stage)

df_mean <- TAX %>% group_by(Compartment_Stage) %>%
    summarise(Mean = mean(Distance))
write.table(df_mean, "Figure_S6c_TAX_mean.txt",
            quote = F, sep = "\t", row.names = F)
lines_mean <- as.numeric(as.character(df_mean$Mean))

p_s6c_2 <- ggplot(TAX, aes(Distance, fill = Compartment_Stage)) +
        geom_density(alpha = 0.7, size = 0.2, color = "wheat4") +
        scale_fill_manual(values = colors) +
        main_theme +
        theme(legend.position = "none", axis.title = element_blank()) +
        geom_vline(xintercept = lines_mean, color = colors,
                   linetype = "longdash")

## put together
p_s6c <- grid.arrange(p_s6c_1, p_s6c_2,
                       p_s6c_1, p_s6c_2,
                       p_s6c_1, p_s6c_2,
                       nrow = 6, ncol = 1)

p <- grid.arrange(p_s6b, p_s6c, nrow = 1, ncol = 2)
ggsave("Figure_S6bc.pdf", p, width = 10, height = 5)
