#!/biodata/dep_psl/grp_rgo/tools/bin/Rscript

library(ggplot2)

c_Man <- data.frame(group = c("NK", "NPK", "CONMIN", "BIODYN"),
                    color = c(c_grey, c_black, c_red, c_green))

c_Gen <- data.frame(group = c("B73", "pht1;6"),
                    color = c(set2[3], set2[8]))

s_Sta <- data.frame(group = c("VG", "RP"),
                    shape = c(1, 16))

s_Man_Sta <- data.frame(group = c("NK_VG", "NK_RP", "NPK_VG", "NPK_RP"),
                        shape = c(2, 17, 1, 16))

get_color_df <- function(x) {
    if (x == "Management") return (c_Man)
    if (x == "Genotype") return (c_Gen)
}

get_shape_df <- function(x) {
    if (x == "Stage") return (s_Sta)
    if (x == "Man_Sta") return (s_Man_Sta)
}

order0 <- c("NK_VG", "NK_RP", "NPK_VG", "NPK_RP")

order <- c("VG_NK_B73",
            "VG_NK_pht1;6",
            "VG_NPK_B73",
            "VG_NPK_pht1;6",
            "RP_NK_B73",
            "RP_NK_pht1;6",
            "RP_NPK_B73",
            "RP_NPK_pht1;6")


order2 <-  c("VG_B73",
            "VG_pht1;6",
            "RP_B73",
            "RP_pht1;6")

box <- function(x, color, shape, group, index) {

    color_df <- get_color_df(color)
    shape_df <- get_shape_df(shape)

    x[[color]] <- factor(x[[color]], levels = color_df$group)
    x[[shape]] <- factor(x[[shape]], levels = shape_df$group)

    p <- ggplot(x, aes_string(x = group, y = index)) +
        geom_boxplot(aes_string(color = color), outlier.shape = NA) +
        geom_jitter(aes_string(color = color, shape = shape),
                    size = 2, alpha = 0.8) +
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
