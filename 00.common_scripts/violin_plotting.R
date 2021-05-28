
violin <- function(x, color, shape, group, index, size = 1.2) {

    color_df <- get_color_df(color)
    shape_df <- get_shape_df(shape)

    p <- ggplot(x, aes_string(x = group, y = index)) +
        geom_violin(aes_string(color = color)) +
        geom_jitter(aes_string(color = color, shape = shape), size = size) +
        scale_colour_manual(values=color_df) +
        scale_shape_manual(values=shape_df) +
        main_theme +
        theme(axis.text.x = element_text(angle = 90, vjust = 0.5,
                                         hjust = 1)) +
        labs(x = "")
 }

