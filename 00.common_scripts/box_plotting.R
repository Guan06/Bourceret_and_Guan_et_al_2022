#!/biodata/dep_psl/grp_rgo/tools/bin/Rscript

library(ggplot2)

get_color_df <- function(x) {
    if (x == "Compartment") return (c_Com)
    if (x == "Management") return (c_Man)
    if (x == "Genotype") return (c_Gen)
    if (x == "Pool") return (c_Pool)
    if (x == "Plot") return (c_Plo)
}

get_shape_df <- function(x) {
    if (x == "Stage") return (s_Sta)
    if (x == "Genotype") return (s_Gen)
    if (x == "Location") return (s_Loc)
    if (x == "Plot") return (s_Plo)
    if (x == "Position") return (s_Pos)
    if (x == "Management") return (s_Man)
    if (x == "Pool") return (s_Pool)
    if (x == "Field") return(s_Fie)
}

box <- function(x, color, shape, group, index) {

    color_df <- get_color_df(color)
    shape_df <- get_shape_df(shape)

    x[[color]] <- factor(x[[color]], levels = color_df$group)
    x[[shape]] <- factor(x[[shape]], levels = shape_df$group)

    p <- ggplot(x, aes_string(x = group, y = index)) +
        geom_boxplot(aes_string(color = color), outlier.shape = NA) +
        geom_jitter(aes_string(color = color, shape = shape), size = 1.2) +
        scale_colour_manual(values=as.character(color_df$color)) +
        scale_shape_manual(values=shape_df$shape) +
        main_theme +
        theme(axis.text.x = element_text(angle = 90, vjust = 0.5,
                                         hjust = 1)) +
        labs(x = "")
 }

box_sig <- function(x, group, compare) {
    lst <- as.character(unique(x[[group]]))
    len <- length(lst)
    sig <- c()
    for (i in 1 : (len - 1 )) {
        for (j in (i + 1) : len) {
            sig_i <- x[x[[group]] == lst[i], ]
            sig_j <- x[x[[group]] == lst[j], ]
            sig_ij <- wilcox.test(sig_i[[compare]], sig_j[[compare]],
                                  alternative = "two.sided")
            this_sig <- c(lst[i], lst[j], sig_ij$p.value)
            sig <- rbind(sig, this_sig)
        }
    }
    colnames(sig) <- c("Group1", "Group2", "Significance")
    rownames(sig) <- c()
    sig <- as.data.frame(sig)
    return(sig)
}
