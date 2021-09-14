#!/biodata/dep_psl/grp_rgo/tools/bin/Rscript

library(randomForest)

get_regressor <- function(asv, response, prefix) {
    ## get the common samples in asv and response(which is design)
    sample_asv <- as.character(colnames(asv))
    sample_res <- as.character(response$Sample_ID)
    inter <- intersect(sample_asv, sample_res)

    asv <- asv[, colnames(asv) %in% inter]
    response <- response[response$Sample_ID %in% inter, ]
    response <- response[match(colnames(asv), response$Sample_ID), ]

    ## filter the ASV with 0 abundance
    asv <- asv[rowSums(asv) > 0, ]
    print(paste0(nrow(asv), " ASVs in ", ncol(asv),
                 " samples will be used as predictor."))

    asv <- t(asv)
    
    num <- ncol(response)
    rf_reg <- list()
    im_mse <- c()

    # skip the first column, which is Sample_ID
    for ( i in 2 : num ) {
        y <- response[, i]
        this_rf <- randomForest(asv, y, ntree = 1000, importance = TRUE)
        rf_reg[[(i - 1)]] <- this_rf

        this_rf_im <- importance(this_rf, scale = TRUE)
        im_mse <- cbind(im_mse, this_rf_im[, 1])
    }
    names(rf_reg) <- colnames(response)[2 : num]
    saveRDS(rf_reg, paste0(prefix, ".rds"))

    colnames(im_mse) <- colnames(response)[2 : num]
    return(im_msell)
}

add_tax <- function(x, tax, level = "Class") {
    id <- rownames(x)

    x_tax <- tax[rownames(tax) %in% id, ]
    x_tax <- x_tax[match(id, rownames(x_tax)), ]

    id_tax <- paste0(id, "_", x_tax[[level]])
    return(id_tax)
}
