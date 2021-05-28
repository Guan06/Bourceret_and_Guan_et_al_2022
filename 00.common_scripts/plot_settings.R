#!/biodata/dep_psl/grp_rgo/tools/bin/Rscript

library(ggplot2, quietly = T, warn.conflicts = F)
library(scales, quietly = T, warn.conflicts = F)
library(grid, quietly = T, warn.conflicts = F)
library(RColorBrewer, quietly = T, warn.conflicts = F)

alpha <- .7
c_yellow <-          rgb(255 / 255, 255 / 255,   0 / 255, alpha)
c_blue <-            rgb(  0 / 255, 000 / 255, 255 / 255, alpha)
c_orange <-          rgb(255 / 255,  69 / 255,   0 / 255, alpha)
c_green <-           rgb(  50/ 255, 220 / 255,  50 / 255, alpha)
c_dark_green <-      rgb( 50 / 255, 200 / 255, 100 / 255, alpha)
c_very_dark_green <- rgb( 50 / 255, 150 / 255, 100 / 255, alpha)
c_sea_green <-       rgb( 46 / 255, 129 / 255,  90 / 255, alpha)
c_black <-           rgb(  0 / 255,   0 / 255,   0 / 255, alpha)
c_grey <-            rgb(180 / 255, 180 / 255,  180 / 255, alpha)
c_dark_brown <-      rgb(101 / 255,  67 / 255,  33 / 255, alpha)
c_red <-             rgb(200 / 255,   0 / 255,   0 / 255, alpha)
c_dark_red <-        rgb(255 / 255, 130 / 255,   0 / 255, alpha)
c_pink <-            rgb(255 / 255,  51 / 255, 204 / 255, alpha)
c_acido <- "#EF5656"
c_bacteroi <- "#47B3DA"
c_firmi <- "#F7A415"
c_alphapro <- "#18E560"
c_betapro <- "#B3EE3A"
c_gammapro <- "#008B00"
c_proteo <- "#2BB065"
c_actino <- "#654321"
c_royalblue1 <- "#4876ff"
c_maroon2 <- "#800000"
c_khaki <- "#C3B091"
c_honeydew <- "#F0FFF0"

set1 <- brewer.pal(n = 11, name = "RdYlGn")
set2 <- brewer.pal(n = 9, name = "Blues")
set3 <- brewer.pal(n = 9, name = "YlOrRd")

## sanzo c312
c_yellow_orange <- "#ff8c00"
c_burnt_sienna <- "#a93400"

## sanzo c249 https://sanzo-wada.dmbk.io/combination/249
c_NK <- c_ecru <- "#c0b490"
c_NPK <- c_hays_russet <- "#681916"
c_CONMIN <- c_olive_ocher <- "#d1bd19"
c_BIODYN <- c_dark_medici_blue <- "#417777"

## set up colors and shapes for different grouping
c_Com <- c("Rhizosphere" = c_dark_red, "Root" = c_very_dark_green,
           "Bulksoil" = c_dark_brown)

c_Man <- c("NK" = c_NK, "NPK" = c_NPK, "CONMIN" = c_CONMIN, "BIODYN" = c_BIODYN)

c_Gen <- c("1_B73" = set2[3], "2_DK105" = set3[3], "3_PH207" = set2[5],
           "4_F2" = set3[4], "5_pht1;6" = set2[8], "Soil" = c_dark_brown)

c_Pool <- c("Dent" = c_blue, "Flint" = c_yellow)

# tints and shades generated from c_CONMIN and c_BIODYN
# https://maketintsandshades.com/#d1bd19,417777
c_Plo <- c("19" = c_NK, "23" = c_NPK,
           "6" = "#e8de8c", "56" = "#d6c430", "96" = "#bcaa17",
           "12" = "#a0bbbb", "50" = "#548585", "90" = "#2e5353")

s_Sta <- c("before_sowing" = 17, "Vegetative" = 1, "Reproductive" = 16)

s_Gen <- c("1_B73" = 15, "2_DK105" = 16, "3_PH207" = 0, "4_F2" = 1,
            "5_pht1;6" = 2, "Soil" = 4)

s_Pool <-c("Dent" = 16, "Flint" = 17)

s_Fie <- c("DEMO" = 17, "DOK" = 2)

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

## theme for the plots
main_theme <- theme(panel.background = element_blank(),
                    panel.grid = element_blank(),
                    plot.title = element_text(hjust = 0.5,
                                              vjust = 0.5,
                                              size = 18),
                    axis.line.x = element_line(color = "black"),
                    axis.line.y = element_line(color = "black"),
                    axis.ticks = element_line(color = "black"),
                    axis.text = element_text(colour = "black", size = 10),
                    legend.position = "right",
                    legend.background = element_blank(),
                    legend.key = element_blank(),
                    text = element_text(family="sans"))

Management_Stage_levels <- c("NK_before_sowing",
                             "NK_Vegetative",
                             "NK_Reproductive",
                             "NPK_before_sowing",
                             "NPK_Vegetative",
                             "NPK_Reproductive",
                             "CONMIN_before_sowing",
                             "CONMIN_Vegetative",
                             "CONMIN_Reproductive",
                             "BIODYN_before_sowing",
                             "BIODYN_Vegetative",
                             "BIODYN_Reproductive")

Plot_Stage_levels <- c("19_before_sowing",
                        "19_Vegetative",
                        "19_Reproductive",
                        "23_before_sowing",
                        "23_Vegetative",
                        "23_Reproductive",
                        "6_before_sowing",
                        "6_Vegetative",
                        "6_Reproductive",
                        "56_before_sowing",
                        "56_Vegetative",
                        "56_Reproductive",
                        "96_before_sowing",
                        "96_Vegetative",
                        "96_Reproductive",
                        "12_before_sowing",
                        "12_Vegetative",
                        "12_Reproductive",
                        "50_before_sowing",
                        "50_Vegetative",
                        "50_Reproductive",
                        "90_before_sowing",
                        "90_Vegetative",
                        "90_Reproductive"
                    )
Plot_Stage_labels <- c("NK_Un", "NK_V", "NK_R",
                       "NPK_Un", "NPK_V", "NPK_R",
                       "CO-1_Un", "CO-1_V", "CO-1_R",
                       "CO-2_Un", "CO-2_V", "CO-2_R",
                       "CO-3_Un", "CO-3_V", "CO-3_R",
                       "BI-1_Un", "BI-1_V", "BI-1_R",
                       "BI-2_Un", "BI-2_V", "BI-2_R",
                       "BI-3_Un", "BI-3_V", "BI-3_R")

Compartment_Stage_levels <- c("Bulksoil_before_sowing",
                              "Bulksoil_Vegetative",
                              "Bulksoil_Reproductive",
                              "Rhizosphere_Vegetative",
                              "Rhizosphere_Reproductive",
                              "Root_Vegetative",
                              "Root_Reproductive")

greens <- brewer.pal(n = 9, name = "Greens")
oranges <- brewer.pal(n = 9, name = "Oranges")
br <- brewer.pal(n = 11, name = "BrBG")

c_Com_Sta <- c("Bulksoil_before_sowing" = br[4],
            "Bulksoil_Vegetative" = br[2],
            "Bulksoil_Reproductive" = br[1],
            "Rhizosphere_Vegetative" = oranges[3],
            "Rhizosphere_Reproductive" = oranges[6],
            "Root_Vegetative" = greens[3],
            "Root_Reproductive" = greens[6])
