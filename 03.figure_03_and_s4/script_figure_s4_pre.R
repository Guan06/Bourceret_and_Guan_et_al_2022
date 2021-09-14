#!/biodata/dep_psl/grp_rgo/tools/bin/Rscript

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

write.table(asv_tax_sum, "Figure_7a_asv_tax_sum.txt",
            quote = F, sep = "\t", row.names = F)
