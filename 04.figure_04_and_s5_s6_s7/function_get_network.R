#!/biodata/dep_psl/grp_rgo/guan/bin/Rscript

library(Hmisc, quietly = T, warn.conflicts = F)
library(igraph, quietly = T, warn.conflicts = F)

###############################################################################
## function modified from http://www.sthda.com/english/wiki/
## correlation-matrix-formatting-and-visualization
flattenCorrMatrix <- function(cormat, pmat) {
    ut <- upper.tri(cormat)
    data.frame(
        node1 = rownames(cormat)[row(cormat)[ut]],
        node2 = rownames(cormat)[col(cormat)[ut]],
        edge  =(cormat)[ut],
        p = pmat[ut]
    )
}

## calculate degree, closeness and betweenness from data frame
get_node_features <- function(x) {
    g <- graph_from_data_frame(x, directed = F)
    dg <- degree(g)
#    cl <- closeness(g)
#    bt <- betweenness(g)

    df <- data.frame(Taxon = names(dg), Degree = dg)
#                     Closeness = cl, Betweenness = bt)
}
###############################################################################
