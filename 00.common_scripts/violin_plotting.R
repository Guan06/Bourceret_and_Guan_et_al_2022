
violin <- function(x, color, shape, group, index) {

    color_df <- get_color_df(color)
    shape_df <- get_shape_df(shape)

    x[[color]] <- factor(x[[color]], levels = color_df$group)
    x[[shape]] <- factor(x[[shape]], levels = shape_df$group)

    p <- ggplot(x, aes_string(x = group, y = index)) +
        geom_violin(aes_string(color = color)) +
        geom_jitter(aes_string(color = color, shape = shape), size = 1.2) +
        scale_colour_manual(values=as.character(color_df$color)) +
        scale_shape_manual(values=shape_df$shape) +
        main_theme +
        theme(axis.text.x = element_text(angle = 90, vjust = 0.5,
                                         hjust = 1)) +
        labs(x = "")
 }

