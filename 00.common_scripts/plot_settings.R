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

## set up colors and shapes for different grouping
c_Com <- data.frame(group = c("Rhizosphere", "Root", "Bulksoil"),
                    color = c(c_dark_red, c_very_dark_green, c_dark_brown))

c_Man <- data.frame(group = c("NK", "NPK", "CONMIN", "BIODYN"),
                    color = c(c_grey, c_black, c_red, c_green))

c_Gen <- data.frame(group = c("1_B73", "2_DK105", "3_PH207", "4_F2",
                            "5_pht1;6", "Soil"),
                    color = c(set2[3], set3[3], set2[5], set3[4], set2[8],
                              c_dark_brown))

c_Pool <- data.frame(group = c("Dent", "Flint"),
                     color = c(c_blue, c_yellow))

c_Plo <- data.frame(group = c("19", "23", "6", "56", "96", "12", "50", "90"),
                    color = c(c_grey, c_black, set1[3], set1[2], set1[1],
                              set1[8], set1[9], set1[11]))

s_Sta <- data.frame(group = c("before_sowing", "Vegetative", "Reproductive"),
                    shape = c(2, 1, 16))

s_Gen <- data.frame(group = c("1_B73", "2_DK105", "3_PH207", "4_F2",
                            "5_pht1;6", "Soil"),
                    shape = c(15, 16, 0, 1, 2, 4))

s_Loc <- data.frame(group = c("A1", "A2", "B1", "B2","B4"),
                    shape = c(1, 2, 17, 16, 15))

s_Plo <- data.frame(group = c("19", "23", "6", "56", "96", "12", "50", "90"),
                    shape = c(2, 1, 15, 16, 17, 15, 16, 17))

s_Pos <- data.frame(group = c("1", "2", "3", "4", "5", "6"),
                    shape = c(0, 1, 2, 15, 16, 17))

s_Man <- data.frame(group = c("NK", "NPK", "CONMIN", "BIODYN"),
                    shape = c(0, 1, 16, 17))

s_Pool <- data.frame(group = c("Dent", "Flint"),
                     shape = c(16, 17))

s_Fie <- data.frame(group = c("DEMO", "DOK"),
                    shape = c(17, 2))

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

Compartment_Stage_levels <- c("Bulksoil_before_sowing",
                              "Bulksoil_Vegetative",
                              "Bulksoil_Reproductive",
                              "Rhizosphere_Vegetative",
                              "Rhizosphere_Reproductive",
                              "Root_Vegetative",
                              "Root_Reproductive")
