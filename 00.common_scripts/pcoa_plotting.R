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

pcoa <- function(x, des, dim, color, shape, size = 1.2, prefix) {
    eig <- x$eig
    eig[eig < 0] <- 0

    if (dim == "12") {
        points <- x$points[, 1 : 2]
    } else if (dim == "34") {
        points <- x$points[, 3 : 4]
    } else if (dim == "13") {
        points <- x$points[, c(1, 3)]
    } else if (dim == "23") {
        points <- x$points[, c(2, 3)]
    }

    eig1 <- 100 * eig[1] / sum(eig)
    eig2 <- 100 * eig[2] / sum(eig)
    eig3 <- 100 * eig[3] / sum(eig)
    eig4 <- 100 * eig[4] / sum(eig)

    points <- as.data.frame(points)
    colnames(points) <- c("x", "y")

    des <- des[match(rownames(points), des[["Sample_ID"]]), ]

    points <- cbind(points, des[[color]])
    colnames(points)[3] <- color
    points <- cbind(points, des[[shape]])
    colnames(points)[4] <- shape

    color_df <- get_color_df(color)
    shape_df <- get_shape_df(shape)

    p <- ggplot(points,
                aes_string(x = "x", y = "y", color = color, shape = shape)) +
         geom_point(size = size, alpha = 0.9) +
         scale_colour_manual(values=color_df) +
         scale_shape_manual(values=shape_df) +
         main_theme 
 #        coord_fixed(ratio = 1)

     if (dim == "12") {
        p <- p + labs (x = paste0("PCo1 (", format(eig1, digits = 4), "%)"),
                       y = paste0("PCo2 (", format(eig2, digits = 4), "%)")) 
#                 ggtitle(paste0(prefix, "_PCo12"))
     } else if (dim == "34") {
        p <- p + labs (x = paste0("PCo3 (", format(eig3, digits = 4), "%)"),
                       y = paste0("PCo4 (", format(eig4, digits = 4), "%)")) 
 #                ggtitle(paste0(prefix, "_PCo34"))
     } else if (dim == "13") {
         p <- p + labs (x = paste0("PCo1 (", format(eig1, digits = 4), "%)"),
                        y = paste0("PCo3 (", format(eig3, digits = 4), "%)")) 
#                 ggtitle(paste0(prefix, "_PCo13"))
     } else if (dim == "23") {
         p <- p + labs (x = paste0("PCo2 (", format(eig2, digits = 4), "%)"),
                        y = paste0("PCo3 (", format(eig3, digits = 4), "%)")) 
 #                ggtitle(paste0(prefix, "_PCo23"))
     }
     return(p)
}
