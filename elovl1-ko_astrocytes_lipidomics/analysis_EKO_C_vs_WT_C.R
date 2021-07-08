library(plyr)
library(tidyverse)

library(readxl)

## For cells
EKO_C_VS_WT_C_cells <- c("EKO-C1-Cells_r1", "EKO-C1-Cells_r2", "EKO-C1-Cells_r3",
                        "EKO-C2-Cells_r1", "EKO-C2-Cells_r2", "EKO-C2-Cells_r3",
                        "EKO-C3-Cells_r1", "EKO-C3-Cells_r2", "EKO-C3-Cells_r3",
                        "WT-C1-Cells_r1", "WT-C1-Cells_r2", "WT-C1-Cells_r3",
                        "WT-C2-Cells_r1", "WT-C2-Cells_r2", "WT-C2-Cells_r3",
                        "WT-C3-Cells_r1", "WT-C3-Cells_r2", "WT-C3-Cells_r3",
                        "blank")

EKO_R_VS_WT_R_cells <- c("EKO-R1-Cells_r1", "EKO-R1-Cells_r2", "EKO-R1-Cells_r3",
                   "EKO-R2-Cells_r1", "EKO-R2-Cells_r2", "EKO-R2-Cells_r3",
                   "EKO-R3-Cells_r1", "EKO-R3-Cells_r2", "EKO-R3-Cells_r3",
                   "WT-R1-Cells_r1", "WT-R1-Cells_r2", "WT-R1-Cells_r3",
                   "WT-R2-Cells_r1", "WT-R2-Cells_r2", "WT-R2-Cells_r3",
                   "WT-R3-Cells_r1", "WT-R3-Cells_r2", "WT-R3-Cells_r3",
                   "blank")

# EKO_R_VS_WT_R_cells <- c("E", "E", "E",
#                          "E", "E", "E",
#                          "E", "E", "E",
#                          "W", "W", "W",
#                          "W", "W", "W",
#                          "W", "W", "W",
#                          "blank")


EKO_C_VS_EKO_R_cells <- c("C1-Cells_r1", "C1-Cells_r2", "C1-Cells_r3",
                   "C2-Cells_r1", "C2-Cells_r2", "C2-Cells_r3",
                   "C3-Cells_r1", "C3-Cells_r2", "C3-Cells_r3",
                   "R1-Cells_r1", "R1-Cells_r2", "R1-Cells_r3",
                   "R2-Cells_r1", "R2-Cells_r2", "R2-Cells_r3",
                   "R3-Cells_r1", "R3-Cells_r2", "R3-Cells_r3",
                   "blank")

WT_C_VS_WT_R_cells <- c("C1-Cells_r1", "C1-Cells_r2", "C1-Cells_r3",
                        "C2-Cells_r1", "C2-Cells_r2", "C2-Cells_r3",
                        "C3-Cells_r1", "C3-Cells_r2", "C3-Cells_r3",
                        "R1-Cells_r1", "R1-Cells_r2", "R1-Cells_r3",
                        "R2-Cells_r1", "R2-Cells_r2", "R2-Cells_r3",
                        "R3-Cells_r1", "R3-Cells_r2", "R3-Cells_r3",
                        "blank")

## For media
EKO_C_VS_WT_C_media <- c("EKO-C1-Media_r1", "EKO-C1-Media_r2", "EKO-C1-Media_r3",
                         "EKO-C2-Media_r1", "EKO-C2-Media_r2", "EKO-C2-Media_r3",
                         "EKO-C3-Media_r1", "EKO-C3-Media_r2", "EKO-C3-Media_r3",
                         "WT-C1-Media_r1", "WT-C1-Media_r2", "WT-C1-Media_r3",
                         "WT-C2-Media_r1", "WT-C2-Media_r2", "WT-C2-Media_r3",
                         "WT-C3-Media_r1", "WT-C3-Media_r2", "WT-C3-Media_r3",
                         "blank")

EKO_R_VS_WT_R_media <- c("EKO-R1-Media_r1", "EKO-R1-Media_r3",
                         "EKO-R2-Media_r1", "EKO-R2-Media_r2", "EKO-R2-Media_r3",
                         "EKO-R3-Media_r1", "EKO-R3-Media_r2", "EKO-R3-Media_r3",
                         "WT-R1-Media_r1", "WT-R1-Media_r3",
                         "WT-R2-Media_r1", "WT-R2-Media_r2", "WT-R2-Media_r3",
                         "WT-R3-Media_r1", "WT-R3-Media_r2", "WT-R3-Media_r3",
                         "blank")


EKO_C_VS_EKO_R_media <- c("C1-Media_r1", "C1-Media_r3",
                          "C2-Media_r1", "C2-Media_r2", "C2-Media_r3",
                          "C3-Media_r1", "C3-Media_r2", "C3-Media_r3",
                          "R1-Media_r1", "R1-Media_r3",
                          "R2-Media_r1", "R2-Media_r2", "R2-Media_r3",
                          "R3-Media_r1", "R3-Media_r2", "R3-Media_r3",
                          "blank")

WT_C_VS_WT_R_media <- c("C1-Media_r1", "C1-Media_r2", "C1-Media_r3",
                        "C2-Media_r1", "C2-Media_r2", "C2-Media_r3",
                        "C3-Media_r1", "C3-Media_r2", "C3-Media_r3",
                        "R1-Media_r1", "R1-Media_r2", "R1-Media_r3",
                        "R2-Media_r1", "R2-Media_r2", "R2-Media_r3",
                        "R3-Media_r1", "R3-Media_r2", "R3-Media_r3",
                        "blank")



lipid_types <- c("FFA" = "FFA", "PC" = "PC")

read_lipids <- function(prefix, marker, exp_names, types, add_blank = NULL) {
  lapply(names(types), function(type){
    
    dir_name = paste0(prefix, type, "/")
    suffix = paste0("_", types[[type]], "_results.xlsx")
    
    all_names <- lapply(exp_names, function(n) {
      file_name = paste0(dir_name, n,   suffix)#"_",marker,
      
      read_xlsx(file_name) %>%
        rename(lipid = `Lipid name`) %>%
        rename(!!n := `Total Intensity`) %>%
        mutate(type = type) %>%
        na.omit
    })
    
    if (!is.null(add_blank)) {
      blank = read_xlsx(paste0(dir_name, add_blank, suffix)) %>%
        rename(lipid = `Lipid name`) %>%
        rename(!!add_blank := `Total Intensity`) %>%
        mutate(type = type) %>%
        na.omit
      
      all_names = append(all_names, list(blank))
    }
    
    Reduce(function(x, y) merge(x, y), all_names)
  }) %>% bind_rows() %>% as_tibble() %>% filter(!grepl("STD_", lipid))
}

read_lipids_EKO <- function(prefix, marker, exp_names, types, add_blank = NULL) {
  lapply(names(types), function(type){
    
    dir_name = paste0(prefix, type, "/EKO/")
    suffix = paste0("_", types[[type]], "_results.xlsx")
    
    all_names <- lapply(exp_names, function(n) {
      file_name = paste0(dir_name, n, suffix)#"_",marker,
      
      read_xlsx(file_name) %>%
        rename(lipid = 'Lipid name') %>%
        rename(!!n := 'Total Intensity') %>%
        mutate(type = type) %>%
        na.omit
    })
    
    if (!is.null(add_blank)) {
      blank = read_xlsx(paste0(dir_name, add_blank, suffix)) %>%
        rename(lipid = `Lipid name`) %>%
        rename(!!add_blank := `Total Intensity`) %>%
        mutate(type = type) %>%
        na.omit
      
      all_names = append(all_names, list(blank))
    }
    
    Reduce(function(x, y) merge(x, y), all_names)
  }) %>% bind_rows() %>% as_tibble() %>% filter(!grepl("STD_", lipid))
}

read_lipids_WT <- function(prefix, marker, exp_names, types, add_blank = NULL) {
  lapply(names(types), function(type){
    
    dir_name = paste0(prefix, type, "/WT/")
    suffix = paste0("_", types[[type]], "_results.xlsx")
    
    all_names <- lapply(exp_names, function(n) {
      file_name = paste0(dir_name, n, suffix)#"_",marker,
      
      read_xlsx(file_name) %>%
        rename(lipid = 'Lipid name') %>%
        rename(!!n := 'Total Intensity') %>%
        mutate(type = type) %>%
        na.omit
    })
    
    if (!is.null(add_blank)) {
      blank = read_xlsx(paste0(dir_name, add_blank, suffix)) %>%
        rename(lipid = 'Lipid name') %>%
        rename(!!add_blank := 'Total Intensity') %>%
        mutate(type = type) %>%
        na.omit
      
      all_names = append(all_names, list(blank))
    }
    
    Reduce(function(x, y) merge(x, y), all_names)
  }) %>% bind_rows() %>% as_tibble() %>% filter(!grepl("STD_", lipid))
}



EKO_C_VS_WT_C_cells_exp = read_lipids(
  "raw_data/", "Cells",
  EKO_C_VS_WT_C_cells, lipid_types
)

EKO_R_VS_WT_R_cells_exp = read_lipids(
  "raw_data/", "Cells",
  EKO_R_VS_WT_R_cells, lipid_types
)


names(EKO_R_VS_WT_R_cells_exp)[names(EKO_R_VS_WT_R_cells_exp) == "EKO-R1-Cells_r1"] <- "E1"
names(EKO_R_VS_WT_R_cells_exp)[names(EKO_R_VS_WT_R_cells_exp) == "EKO-R1-Cells_r2"] <- "E2"
names(EKO_R_VS_WT_R_cells_exp)[names(EKO_R_VS_WT_R_cells_exp) == "EKO-R1-Cells_r3"] <- "E3"
names(EKO_R_VS_WT_R_cells_exp)[names(EKO_R_VS_WT_R_cells_exp) == "EKO-R2-Cells_r1"] <- "E4"
names(EKO_R_VS_WT_R_cells_exp)[names(EKO_R_VS_WT_R_cells_exp) == "EKO-R2-Cells_r2"] <- "E5"
names(EKO_R_VS_WT_R_cells_exp)[names(EKO_R_VS_WT_R_cells_exp) == "EKO-R2-Cells_r3"] <- "E6"
names(EKO_R_VS_WT_R_cells_exp)[names(EKO_R_VS_WT_R_cells_exp) == "EKO-R3-Cells_r1"] <- "E7"
names(EKO_R_VS_WT_R_cells_exp)[names(EKO_R_VS_WT_R_cells_exp) == "EKO-R3-Cells_r2"] <- "E8"
names(EKO_R_VS_WT_R_cells_exp)[names(EKO_R_VS_WT_R_cells_exp) == "EKO-R3-Cells_r3"] <- "E9"

names(EKO_R_VS_WT_R_cells_exp)[names(EKO_R_VS_WT_R_cells_exp) == "WT-R1-Cells_r1"] <- "W1"
names(EKO_R_VS_WT_R_cells_exp)[names(EKO_R_VS_WT_R_cells_exp) == "WT-R1-Cells_r2"] <- "W2"
names(EKO_R_VS_WT_R_cells_exp)[names(EKO_R_VS_WT_R_cells_exp) == "WT-R1-Cells_r3"] <- "W3"
names(EKO_R_VS_WT_R_cells_exp)[names(EKO_R_VS_WT_R_cells_exp) == "WT-R2-Cells_r1"] <- "W4"
names(EKO_R_VS_WT_R_cells_exp)[names(EKO_R_VS_WT_R_cells_exp) == "WT-R2-Cells_r2"] <- "W5"
names(EKO_R_VS_WT_R_cells_exp)[names(EKO_R_VS_WT_R_cells_exp) == "WT-R2-Cells_r3"] <- "W6"
names(EKO_R_VS_WT_R_cells_exp)[names(EKO_R_VS_WT_R_cells_exp) == "WT-R3-Cells_r1"] <- "W7"
names(EKO_R_VS_WT_R_cells_exp)[names(EKO_R_VS_WT_R_cells_exp) == "WT-R3-Cells_r2"] <- "W8"
names(EKO_R_VS_WT_R_cells_exp)[names(EKO_R_VS_WT_R_cells_exp) == "WT-R3-Cells_r3"] <- "W9"

# EKO_R_VS_WT_R_cells_exp$e1=rowMeans(EKO_R_VS_WT_R_cells_exp[,c("E1", "E2", "E3")], na.rm=TRUE)
# EKO_R_VS_WT_R_cells_exp$e2=rowMeans(EKO_R_VS_WT_R_cells_exp[,c("E4", "E5", "E6")], na.rm=TRUE)
# EKO_R_VS_WT_R_cells_exp$e3=rowMeans(EKO_R_VS_WT_R_cells_exp[,c("E7", "E8", "E9")], na.rm=TRUE)
# 
# EKO_R_VS_WT_R_cells_exp$w1=rowMeans(EKO_R_VS_WT_R_cells_exp[,c("W1", "W2", "W3")], na.rm=TRUE)
# EKO_R_VS_WT_R_cells_exp$w2=rowMeans(EKO_R_VS_WT_R_cells_exp[,c("W4", "W5", "W6")], na.rm=TRUE)
# EKO_R_VS_WT_R_cells_exp$w3=rowMeans(EKO_R_VS_WT_R_cells_exp[,c("W7", "W8", "W9")], na.rm=TRUE)
# 
# EKO_R_VS_WT_R_cells_exp <- select(EKO_R_VS_WT_R_cells_exp, -E1, -E2, -E3, -E4, -E5, -E6, -E7, -E8, -E9,
#                                    -W1, -W2, -W3, -W4, -W5, -W6, -W7, -W8, -W9)

EKO_C_VS_WT_C_media_exp = read_lipids(
  "raw_data/", "Media",
  EKO_C_VS_WT_C_media, lipid_types
)

EKO_R_VS_WT_R_media_exp = read_lipids(
  "raw_data/", "Media",
  EKO_R_VS_WT_R_media, lipid_types
)

names(EKO_R_VS_WT_R_media_exp)[names(EKO_R_VS_WT_R_media_exp) == "EKO-R1-Media_r1"] <- "E1"
names(EKO_R_VS_WT_R_media_exp)[names(EKO_R_VS_WT_R_media_exp) == "EKO-R1-Media_r2"] <- "E2"
names(EKO_R_VS_WT_R_media_exp)[names(EKO_R_VS_WT_R_media_exp) == "EKO-R1-Media_r3"] <- "E3"
names(EKO_R_VS_WT_R_media_exp)[names(EKO_R_VS_WT_R_media_exp) == "EKO-R2-Media_r1"] <- "E4"
names(EKO_R_VS_WT_R_media_exp)[names(EKO_R_VS_WT_R_media_exp) == "EKO-R2-Media_r2"] <- "E5"
names(EKO_R_VS_WT_R_media_exp)[names(EKO_R_VS_WT_R_media_exp) == "EKO-R2-Media_r3"] <- "E6"
names(EKO_R_VS_WT_R_media_exp)[names(EKO_R_VS_WT_R_media_exp) == "EKO-R3-Media_r1"] <- "E7"
names(EKO_R_VS_WT_R_media_exp)[names(EKO_R_VS_WT_R_media_exp) == "EKO-R3-Media_r2"] <- "E8"
names(EKO_R_VS_WT_R_media_exp)[names(EKO_R_VS_WT_R_media_exp) == "EKO-R3-Media_r3"] <- "E9"

names(EKO_R_VS_WT_R_media_exp)[names(EKO_R_VS_WT_R_media_exp) == "WT-R1-Media_r1"] <- "W1"
names(EKO_R_VS_WT_R_media_exp)[names(EKO_R_VS_WT_R_media_exp) == "WT-R1-Media_r2"] <- "W2"
names(EKO_R_VS_WT_R_media_exp)[names(EKO_R_VS_WT_R_media_exp) == "WT-R1-Media_r3"] <- "W3"
names(EKO_R_VS_WT_R_media_exp)[names(EKO_R_VS_WT_R_media_exp) == "WT-R2-Media_r1"] <- "W4"
names(EKO_R_VS_WT_R_media_exp)[names(EKO_R_VS_WT_R_media_exp) == "WT-R2-Media_r2"] <- "W5"
names(EKO_R_VS_WT_R_media_exp)[names(EKO_R_VS_WT_R_media_exp) == "WT-R2-Media_r3"] <- "W6"
names(EKO_R_VS_WT_R_media_exp)[names(EKO_R_VS_WT_R_media_exp) == "WT-R3-Media_r1"] <- "W7"
names(EKO_R_VS_WT_R_media_exp)[names(EKO_R_VS_WT_R_media_exp) == "WT-R3-Media_r2"] <- "W8"
names(EKO_R_VS_WT_R_media_exp)[names(EKO_R_VS_WT_R_media_exp) == "WT-R3-Media_r3"] <- "W9"


EKO_C_VS_EKO_R_cells_exp = read_lipids_EKO(
  "raw_data/", "Media",
  EKO_C_VS_EKO_R_cells, lipid_types
)

EKO_C_VS_EKO_R_media_exp = read_lipids_EKO(
  "raw_data/", "Media",
  EKO_C_VS_EKO_R_media, lipid_types
)

WT_C_VS_WT_R_cells_exp = read_lipids_WT(
  "raw_data/", "Media",
  WT_C_VS_WT_R_cells, lipid_types
)

WT_C_VS_WT_R_media_exp = read_lipids_WT(
  "raw_data/", "Media",
  WT_C_VS_WT_R_media, lipid_types
)


library(edgeR)

# Groups for experiment 1
gr_EKO_C_VS_WT_C_cells_expr = c("EKO_C", "EKO_C","EKO_C","EKO_C","EKO_C","EKO_C","EKO_C","EKO_C","EKO_C",
              "WT_C", "WT_C", "WT_C","WT_C","WT_C","WT_C","WT_C","WT_C","WT_C",
              "Blank") %>%
  factor(levels = c("Blank", "EKO_C", "WT_C"))

design_EKO_C_VS_WT_C_cells = model.matrix(~gr_EKO_C_VS_WT_C_cells_expr)

contrasts_EKO_C_VS_WT_C_cells = makeContrasts(
  H = gr_EKO_C_VS_WT_C_cells_exprEKO_C - gr_EKO_C_VS_WT_C_cells_exprWT_C,
  levels = design_EKO_C_VS_WT_C_cells
)

# Groups for experiment 2
gr_EKO_R_VS_WT_R_cells_expr = c("EKO_R", "EKO_R","EKO_R","EKO_R","EKO_R","EKO_R","EKO_R","EKO_R","EKO_R",
                                "WT_R", "WT_R", "WT_R","WT_R","WT_R","WT_R","WT_R","WT_R","WT_R",
                                "Blank") %>%
  factor(levels = c("Blank", "EKO_R", "WT_R"))

design_EKO_R_VS_WT_R_cells = model.matrix(~gr_EKO_R_VS_WT_R_cells_expr)

contrasts_EKO_R_VS_WT_R_cells = makeContrasts(
  H = gr_EKO_R_VS_WT_R_cells_exprEKO_R - gr_EKO_R_VS_WT_R_cells_exprWT_R,
  levels = design_EKO_R_VS_WT_R_cells
)

# # Groups for experiment 2 means
# gr_EKO_R_VS_WT_R_cells_expr = c("E", "E", "E",
#                                 "W", "W", "W",
#                                 "Blank") %>%
#   factor(levels = c("Blank", "E", "W"))
# 
# design_EKO_R_VS_WT_R_cells = model.matrix(~gr_EKO_R_VS_WT_R_cells_expr)
# 
# contrasts_EKO_R_VS_WT_R_cells = makeContrasts(
#   H = gr_EKO_R_VS_WT_R_cells_exprE - gr_EKO_R_VS_WT_R_cells_exprW,
#   levels = design_EKO_R_VS_WT_R_cells
# )

# Groups for experiment 3
gr_EKO_C_VS_WT_C_media_expr = c("EKO_C", "EKO_C","EKO_C","EKO_C","EKO_C","EKO_C","EKO_C","EKO_C","EKO_C",
                                "WT_C", "WT_C", "WT_C","WT_C","WT_C","WT_C","WT_C","WT_C","WT_C",
                                "Blank") %>%
  factor(levels = c("Blank", "EKO_C", "WT_C"))

design_EKO_C_VS_WT_C_media = model.matrix(~gr_EKO_C_VS_WT_C_media_expr)

contrasts_EKO_C_VS_WT_C_media = makeContrasts(
  H = gr_EKO_C_VS_WT_C_media_exprEKO_C - gr_EKO_C_VS_WT_C_media_exprWT_C,
  levels = design_EKO_C_VS_WT_C_media
)

# Groups for experiment 4


gr_EKO_R_VS_WT_R_media_expr = c("EKO_R", "EKO_R","EKO_R","EKO_R","EKO_R","EKO_R","EKO_R","EKO_R",
                                "WT_R", "WT_R", "WT_R","WT_R","WT_R","WT_R","WT_R","WT_R",
                                "Blank") %>%
  factor(levels = c("Blank", "EKO_R", "WT_R"))

design_EKO_R_VS_WT_R_media = model.matrix(~gr_EKO_R_VS_WT_R_media_expr)

contrasts_EKO_R_VS_WT_R_media = makeContrasts(
  H = gr_EKO_R_VS_WT_R_media_exprEKO_R - gr_EKO_R_VS_WT_R_media_exprWT_R,
  levels = design_EKO_R_VS_WT_R_media
)


# Groups for experiment 5


gr_EKO_C_VS_EKO_R_cells_expr = c("C", "C","C","C","C","C","C","C","C",
                                "R", "R", "R","R","R","R","R","R","R",
                                "Blank") %>%
  factor(levels = c("Blank", "C", "R"))

design_EKO_C_VS_EKO_R_cells = model.matrix(~gr_EKO_C_VS_EKO_R_cells_expr)

contrasts_EKO_C_VS_EKO_R_cells = makeContrasts(
  H = gr_EKO_C_VS_EKO_R_cells_exprC - gr_EKO_C_VS_EKO_R_cells_exprR,
  levels = design_EKO_C_VS_EKO_R_cells
)

# Groups for experiment 6

gr_WT_C_VS_WT_R_cells_expr = c("C", "C","C","C","C","C","C","C","C",
                                 "R", "R", "R","R","R","R","R","R","R",
                                 "Blank") %>%
  factor(levels = c("Blank", "C", "R"))

design_WT_C_VS_WT_R_cells = model.matrix(~gr_WT_C_VS_WT_R_cells_expr)

contrasts_WT_C_VS_WT_R_cells = makeContrasts(
  H = gr_WT_C_VS_WT_R_cells_exprC - gr_WT_C_VS_WT_R_cells_exprR,
  levels = design_WT_C_VS_WT_R_cells
)

# Groups for experiment 7

gr_EKO_C_VS_EKO_R_media_expr = c("C", "C","C","C","C","C","C","C",
                               "R", "R", "R","R","R","R","R","R",
                               "Blank") %>%
  factor(levels = c("Blank", "C", "R"))

design_EKO_C_VS_EKO_R_media = model.matrix(~gr_EKO_C_VS_EKO_R_media_expr)

contrasts_EKO_C_VS_EKO_R_media = makeContrasts(
  H = gr_EKO_C_VS_EKO_R_media_exprC - gr_EKO_C_VS_EKO_R_media_exprR,
  levels = design_EKO_C_VS_EKO_R_media
)



# Groups for experiment 8

gr_WT_C_VS_WT_R_media_expr = c("C", "C","C","C","C","C","C","C","C",
                                 "R", "R", "R","R","R","R","R","R","R",
                                 "Blank") %>%
  factor(levels = c("Blank", "C", "R"))

design_WT_C_VS_WT_R_media = model.matrix(~gr_WT_C_VS_WT_R_media_expr)

contrasts_WT_C_VS_WT_R_media = makeContrasts(
  H = gr_WT_C_VS_WT_R_media_exprC - gr_WT_C_VS_WT_R_media_exprR,
  levels = design_WT_C_VS_WT_R_media
)


#######################################

perform_analysis_raw <- function(counts, design_mat, gr) {
  
  data.edgeR <- DGEList(counts = counts %>%
                          na.omit %>%
                          mutate(lipid = make.unique(lipid)) %>%
                          select(-Transition, -type) %>%
                          column_to_rownames("lipid"),
                        group = gr
  )
  
  data.edgeR <- calcNormFactors(data.edgeR, method="TMM")
  data.edgeR <- estimateCommonDisp(data.edgeR, design=design_mat)
  data.edgeR
}

calculate_significance <- function(dge, contrast) {
  dge %>%
    glmFit() %>%
    glmLRT(contrast = contrast)
}

EKO_C_VS_WT_C_cells_expr_helper <- function(df) {
  df %>%
    perform_analysis_raw(design_EKO_C_VS_WT_C_cells, gr_EKO_C_VS_WT_C_cells_expr) %>%
    calculate_significance(contrasts_EKO_C_VS_WT_C_cells) %>%
    topTags(1500000) %>%
    as.data.frame() %>%
    rownames_to_column("lipid") %>%
    as_tibble()
  
}

EKO_R_VS_WT_R_cells_expr_helper <- function(df) {
  df %>%
    perform_analysis_raw(design_EKO_R_VS_WT_R_cells, gr_EKO_R_VS_WT_R_cells_expr) %>%
    calculate_significance(contrasts_EKO_R_VS_WT_R_cells) %>%
    topTags(1500000) %>%
    as.data.frame() %>%
    rownames_to_column("lipid") %>%
    as_tibble()
  
}


EKO_C_VS_WT_C_media_expr_helper <- function(df) {
  df %>%
    perform_analysis_raw(design_EKO_C_VS_WT_C_media, gr_EKO_C_VS_WT_C_media_expr) %>%
    calculate_significance(contrasts_EKO_C_VS_WT_C_media) %>%
    topTags(1500000) %>%
    as.data.frame() %>%
    rownames_to_column("lipid") %>%
    as_tibble()
}


EKO_R_VS_WT_R_media_expr_helper <- function(df) {
  df %>%
    perform_analysis_raw(design_EKO_R_VS_WT_R_media, gr_EKO_R_VS_WT_R_media_expr) %>%
    calculate_significance(contrasts_EKO_R_VS_WT_R_media) %>%
    topTags(1500000) %>%
    as.data.frame() %>%
    rownames_to_column("lipid") %>%
    as_tibble()
}


EKO_C_VS_EKO_R_cells_expr_helper <- function(df) {
  df %>%
    perform_analysis_raw(design_EKO_C_VS_EKO_R_cells, gr_EKO_C_VS_EKO_R_cells_expr) %>%
    calculate_significance(contrasts_EKO_C_VS_EKO_R_cells) %>%
    topTags(1500000) %>%
    as.data.frame() %>%
    rownames_to_column("lipid") %>%
    as_tibble()
}



WT_C_VS_WT_R_cells_expr_helper <- function(df) {
  df %>%
    perform_analysis_raw(design_WT_C_VS_WT_R_cells, gr_WT_C_VS_WT_R_cells_expr) %>%
    calculate_significance(contrasts_WT_C_VS_WT_R_cells) %>%
    topTags(1500000) %>%
    as.data.frame() %>%
    rownames_to_column("lipid") %>%
    as_tibble()
}


EKO_C_VS_EKO_R_media_expr_helper <- function(df) {
  df %>%
    perform_analysis_raw(design_EKO_C_VS_EKO_R_media, gr_EKO_C_VS_EKO_R_media_expr) %>%
    calculate_significance(contrasts_EKO_C_VS_EKO_R_media) %>%
    topTags(1500000) %>%
    as.data.frame() %>%
    rownames_to_column("lipid") %>%
    as_tibble()
}


WT_C_VS_WT_R_media_expr_helper <- function(df) {
  df %>%
    perform_analysis_raw(design_WT_C_VS_WT_R_media, gr_WT_C_VS_WT_R_media_expr) %>%
    calculate_significance(contrasts_WT_C_VS_WT_R_media) %>%
    topTags(1500000) %>%
    as.data.frame() %>%
    rownames_to_column("lipid") %>%
    as_tibble()
}



EKO_C_VS_WT_C_cells_exp_tbl <-
  EKO_C_VS_WT_C_cells_exp %>%
  EKO_C_VS_WT_C_cells_expr_helper

EKO_R_VS_WT_R_cells_exp_tbl <-
  EKO_R_VS_WT_R_cells_exp %>%
  EKO_R_VS_WT_R_cells_expr_helper

EKO_C_VS_WT_C_media_exp_tbl <-
  EKO_C_VS_WT_C_media_exp %>%
  EKO_C_VS_WT_C_media_expr_helper

EKO_R_VS_WT_R_media_exp_tbl <-
  EKO_R_VS_WT_R_media_exp %>%
  EKO_R_VS_WT_R_media_expr_helper


EKO_C_VS_EKO_R_cells_exp_tbl <-
  EKO_C_VS_EKO_R_cells_exp %>%
  EKO_C_VS_EKO_R_cells_expr_helper


WT_C_VS_WT_R_cells_exp_tbl <-
  WT_C_VS_WT_R_cells_exp %>%
  WT_C_VS_WT_R_cells_expr_helper

EKO_C_VS_EKO_R_media_exp_tbl <-
  EKO_C_VS_EKO_R_media_exp %>%
  EKO_C_VS_EKO_R_media_expr_helper


WT_C_VS_WT_R_media_exp_tbl <-
  WT_C_VS_WT_R_media_exp %>%
  WT_C_VS_WT_R_media_expr_helper


EKO_C_VS_WT_C_cells_exp_tbl

EKO_R_VS_WT_R_cells_exp_tbl

EKO_C_VS_WT_C_media_exp_tbl

EKO_R_VS_WT_R_media_exp_tbl

EKO_C_VS_EKO_R_cells_exp_tbl

WT_C_VS_WT_R_cells_exp_tbl

EKO_C_VS_EKO_R_media_exp_tbl

WT_C_VS_WT_R_media_exp_tbl


make_volcano_plot <- function(df, title) {
  df %>%
    mutate(sig = factor(FDR < 0.10)) %>%
    ggplot(aes(logFC, -log10(FDR), color = sig)) +
    geom_point() +
    scale_color_manual(values = c("FALSE" = "black", "TRUE" = "red")) +
    guides(color = F) +
    ggtitle(title)
}

EKO_C_VS_WT_C_cells_exp_tbl %>% make_volcano_plot("EKO Control cells vs WT Control cells")

EKO_R_VS_WT_R_cells_exp_tbl %>% make_volcano_plot("EKO Reactive cells vs WT Reactive cells")

EKO_C_VS_WT_C_media_exp_tbl %>% make_volcano_plot("EKO Control Media vs WT Control Media")

EKO_R_VS_WT_R_media_exp_tbl %>% make_volcano_plot("EKO Reactive Media vs WT Reactive Media")

EKO_C_VS_EKO_R_cells_exp_tbl %>% make_volcano_plot("EKO Control cells vs EKO Reactive cells")

WT_C_VS_WT_R_cells_exp_tbl %>% make_volcano_plot("WT Control cells vs WT Reactive cells")

EKO_C_VS_EKO_R_media_exp_tbl %>% make_volcano_plot("EKO Control Media vs EKO Reactive Media")

WT_C_VS_WT_R_media_exp_tbl %>% make_volcano_plot("WT Control Media vs WT Reactive Media")


write_summary_and_results <- function(tbl, df, name) {
  
  tbl %>%
    merge(df) %>%
    as_tibble() %>%
    arrange(FDR) -> results
  
  write_csv(results, paste0("results/", name, "_full.csv"))
  
  results %>%
    group_by(type) %>%
    summarise(ab = sum(logFC < 0 & FDR < 0.1),
              ve = sum(logFC > 0 & FDR < 0.1)) %>%
    write_csv(paste0("results/", name, "_summary.csv"))
}

dir.create("results", F)

EKO_C_VS_WT_C_cells_exp_tbl %>% write_summary_and_results(EKO_C_VS_WT_C_cells_exp, "EKO_C_VS_WT_C_cells")
EKO_R_VS_WT_R_cells_exp_tbl %>% write_summary_and_results(EKO_R_VS_WT_R_cells_exp, "EKO_R_VS_WT_R_cells")
EKO_C_VS_WT_C_media_exp_tbl %>% write_summary_and_results(EKO_C_VS_WT_C_media_exp, "EKO_C_VS_WT_C_media")
EKO_R_VS_WT_R_media_exp_tbl %>% write_summary_and_results(EKO_R_VS_WT_R_media_exp, "EKO_R_VS_WT_R_media")
EKO_C_VS_EKO_R_cells_exp_tbl %>% write_summary_and_results(EKO_C_VS_EKO_R_cells_exp, "EKO_C_VS_EKO_R_cells")
WT_C_VS_WT_R_cells_exp_tbl %>% write_summary_and_results(WT_C_VS_WT_R_cells_exp, "WT_C_VS_WT_R_cells")
EKO_C_VS_EKO_R_media_exp_tbl %>% write_summary_and_results(EKO_C_VS_EKO_R_media_exp, "EKO_C_VS_EKO_R_media")
WT_C_VS_WT_R_media_exp_tbl %>% write_summary_and_results(WT_C_VS_WT_R_media_exp, "WT_C_VS_WT_R_media")


source("ggbiplot.R")

get_DE_lipids <- function(counts, design_mat, gr, contrasts, p.value = 0.1) {
  dls <-
    counts %>%
    perform_analysis_raw(design_mat, gr) %>%
    calculate_significance(contrasts) %>%
    decideTestsDGE(p.value = p.value)
  
  rownames(dls)[dls %>% as.logical()]
}

make_pca_plot <- function(tp, design_mat, gr, contrasts,
                          title = "PCA plot",
                          ellipse = T, var.axes = F,
                          labels = T) {
  
  tp_edger <-
    tp %>%
    get_DE_lipids(design_mat, gr, contrasts)
  
  if(length(tp_edger) == 0) {
    cat("No significant lipids for ", title)
    return()
  }
  
  if(length(tp_edger) == 1) {
    cat("Single significant lipids for ", title, " is ", tp_edger[1])
    return()
  }
  
  tp %>%
    na.omit %>%
    mutate(lipid = make.unique(lipid)) %>%
    filter(lipid %in% tp_edger) %>%
    select(-Transition, -type, -blank) %>%
    column_to_rownames("lipid") %>%
    as.matrix() %>%
    t %>%
    prcomp(center = T, scale = T) ->
    prcomp_data
  
  groups = NULL
  
  tp %>%
    select(-Transition, -type, -blank, -lipid) %>%
    colnames() ->
    labels.tmp
  
  groups = substr(labels.tmp, 1, 1)
  
  if (!is.null(labels)) {
    labels = labels.tmp
  }
  
  prcomp_data %>%
    ggbiplot(ellipse = ellipse,
             labels = labels,
             groups = groups,
             var.axes = var.axes
    ) + xlim(-2.5,2.5) + ylim(-2.5,2.5)+
    ggtitle(title) +
    cowplot::theme_cowplot()
}

EKO_C_VS_WT_C_cells_exp %>%
  make_pca_plot(design_EKO_C_VS_WT_C_cells, gr_EKO_C_VS_WT_C_cells_expr, contrasts_EKO_C_VS_WT_C_cells, "EKO Control Cells vs WT Control Cells")


EKO_R_VS_WT_R_cells_exp %>%
  make_pca_plot(design_EKO_R_VS_WT_R_cells, gr_EKO_R_VS_WT_R_cells_expr, contrasts_EKO_R_VS_WT_R_cells, "EKO Reactive Cells vs WT Reactive Cells")

EKO_C_VS_WT_C_media_exp %>%
  make_pca_plot(design_EKO_C_VS_WT_C_media, gr_EKO_C_VS_WT_C_media_expr, contrasts_EKO_C_VS_WT_C_media, "EKO Control Media vs WT Control Media")


EKO_R_VS_WT_R_media_exp %>%
  make_pca_plot(design_EKO_R_VS_WT_R_media, gr_EKO_R_VS_WT_R_media_expr, contrasts_EKO_R_VS_WT_R_media, "EKO Reactive Media vs WT Reactive Media")

EKO_C_VS_EKO_R_cells_exp %>%
  make_pca_plot(design_EKO_C_VS_EKO_R_cells, gr_EKO_C_VS_EKO_R_cells_expr, contrasts_EKO_C_VS_EKO_R_cells, "EKO Control Cells vs EKO Reactive Cells")

WT_C_VS_WT_R_cells_exp %>%
  make_pca_plot(design_WT_C_VS_WT_R_cells, gr_WT_C_VS_WT_R_cells_expr, contrasts_WT_C_VS_WT_R_cells, "WT Control Cells vs WT Reactive Cells")

EKO_C_VS_EKO_R_media_exp %>%
  make_pca_plot(design_EKO_C_VS_EKO_R_media, gr_EKO_C_VS_EKO_R_media_expr, contrasts_EKO_C_VS_EKO_R_media, "EKO Control Media vs EKO Reactive Media")


WT_C_VS_WT_R_media_exp %>%
  make_pca_plot(design_WT_C_VS_WT_R_media, gr_WT_C_VS_WT_R_media_expr, contrasts_WT_C_VS_WT_R_media, "WT Control Media vs WT Reactive Media")

