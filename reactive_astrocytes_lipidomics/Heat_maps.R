library(plyr)
library(tidyverse)
library(readxl)

experiment_1_names <- c("C1", "C2", "C3", "C4", "C5",
                        "R1", "R2", "R3", "R4", "R5",
                        "blank")

experiment_2_names <- c("Control1", "Control2", "Control3", "Control4",
                        "RAC1", "RAC2", "RAC3", "RAC4",
                        "blank")

lipid_types <- c("AC" = "AC", "CE" = "CE", "cer" = "cer", "FFA" = "FFA",
                 "SM" = "LysoPC_PC_SM", "PE" = "Lyso_PE", "PG" = "Lyso_PG",
                 "PI" = "Lyso_PI", "PS" = "Lyso_PS",
                 "TAG1" = "TAG1", "TAG2" = "TAG2")

metob_types <- c("NEG1" = "NEG1", "NEG2" = "NEG2", "POS1" = "POS1", "POS2" = "POS2")

read_lipids <- function(prefix, marker, exp_names, types, add_blank = NULL) {
  lapply(names(types), function(type){
    
    dir_name = paste0(prefix, "/", type, "/")
    suffix = paste0("_", types[[type]], "_results.xlsx")
    
    all_names <- lapply(exp_names, function(n) {
      file_name = paste0(dir_name, n, "_", marker, suffix)
      
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

cells_lipid_expr_1 = read_lipids(
  "raw_data/cells_lipids/", "cells",
  experiment_1_names, lipid_types
)

media_lipid_expr_1 = read_lipids(
  "raw_data/media_lipids/", "media", experiment_1_names, lipid_types,
  "blank_solvent"
)

cells_metab_expr_1 = read_lipids(
  "raw_data/cells_and_media_metabolites/", "cells",
  experiment_1_names, metob_types
)

media_metab_expr_1 = read_lipids(
  "raw_data/cells_and_media_metabolites/", "media",
  experiment_1_names, metob_types,
  "solvent_blank_media"
)

cells_lipid_expr_2 = read_lipids(
  "raw_data/cells_lipids/", "cells",
  experiment_2_names, lipid_types
)

media_lipid_expr_2 = read_lipids(
  "raw_data/media_lipids/", "media",
  experiment_2_names, lipid_types, "blank_solvent"
)

cells_metab_expr_2 = read_lipids(
  "raw_data/cells_and_media_metabolites/", "cells",
  experiment_2_names, metob_types
)

media_metab_expr_2 = read_lipids(
  "raw_data/cells_and_media_metabolites/", "media",
  experiment_2_names, metob_types,
  "solvent_blank_media"
)


library(edgeR)

# Groups for experiment 1
gr_expr_1 = c("C", "C", "C", "C", "C",
              "R", "R", "R", "R", "R",
              "Blank") %>%
  factor(levels = c("Blank", "C", "R"))

design_expr_1 = model.matrix(~gr_expr_1)

contrasts_expr_1 = makeContrasts(
  H = gr_expr_1R - gr_expr_1C,
  levels = design_expr_1
)

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

experiment_1_helper <- function(df) {
  df %>%
    perform_analysis_raw(design_expr_1, gr_expr_1) %>%
    calculate_significance(contrasts_expr_1) %>%
    topTags(1500000) %>%
    as.data.frame() %>%
    rownames_to_column("lipid") %>%
    as_tibble()
  
}

cl_e1_tbl <-
  cells_lipid_expr_1 %>%
  experiment_1_helper

ml_e1_tbl <-
  media_lipid_expr_1 %>%
  select(-blank_solvent) %>%
  experiment_1_helper

cm_e1_tbl <-
  cells_metab_expr_1 %>%
  experiment_1_helper

mm_e1_tbl <-
  media_metab_expr_1 %>%
  select(-solvent_blank_media) %>%
  experiment_1_helper


cl_e1_tbl
ml_e1_tbl
cm_e1_tbl
mm_e1_tbl

#######################
# Heat Maps

get_DE_lipids <- function(counts, design_mat, gr, contrasts, p.value = 0.1) {
  dls <-
    counts %>%
    perform_analysis_raw(design_mat, gr) %>%
    calculate_significance(contrasts) %>%
    decideTestsDGE(p.value = p.value)
  
  rownames(dls)[dls %>% as.logical()]
}

library(RColorBrewer)

make_heatmap_jf_e1 <- function(tp, design_mat, gr, contrasts, title = "Heat-map", FDR = 0.10, breaks, list1, list2) {
  
  DElist <-
    tp %>%
    get_DE_lipids(design_mat, gr, contrasts, FDR)
  
  if(length(DElist) == 0) {
    cat("No significant lipids for ", title)
    return()
  }
  
  if(length(DElist) == 1) {
    cat("Single significant lipids for ", title, " is ", DElist[1])
    return()
  }
  
  tp %>%
    mutate(lipid = make.unique(lipid)) %>%
    filter(lipid %in% DElist) %>%
    select(-Transition, -type) %>%
    rowwise() %>%
    mutate(mean = mean(c(C1, C2, C3, C4, C5, R1, R2, R3, R4, R5))) %>%
    mutate_if(is.numeric, log2) %>%
    mutate_if(is.numeric, list(~ . - mean)) %>%
    select(-blank, -mean) ->
    vals
  
  vals %>%
    column_to_rownames("lipid") %>%
    as.matrix() %>%
    pheatmap::pheatmap(color = colorRampPalette(rev(brewer.pal(n = 7, name =
                                                                 "RdBu")))(100), 
                       main = title, 
                       cluster_rows = FALSE,
                       cluster_cols = FALSE, breaks = breaks, legend_breaks = list1, legend_labels = list2,
                       fontsize = 10, fontsize_row = 6)
}

cells_lipid_expr_1 %>%
  make_heatmap_jf_e1(design_expr_1, gr_expr_1, contrasts_expr_1, "Cells Lipids in Experiment 1", 0.001, 
                     seq(-2, 2, length.out = 100), -2:2, -2:2)

media_lipid_expr_1 %>% select(-blank_solvent) %>%
  make_heatmap_jf_e1(design_expr_1, gr_expr_1, contrasts_expr_1, "Media Lipids in Experiment 1", 0.001,
                     seq(-1, 1, length.out = 100), -1:1, -1:1)

cells_metab_expr_1 %>%
  make_heatmap_jf_e1(design_expr_1, gr_expr_1, contrasts_expr_1, "Cells Metabolites in Experiment 1", 0.001,
                     seq(-1.5, 1.5, length.out = 100), -1.5:1.5, -1.5:1.5)

media_metab_expr_1 %>% select(-solvent_blank_media) %>%
  make_heatmap_jf_e1(design_expr_1, gr_expr_1, contrasts_expr_1, "Media Metabolites in Experiment 1", 0.001,
                     seq(-1.3, 1.3, length.out = 100), seq(-1.3, 1.3, length.out = 3), seq(-1.3, 1.3, length.out = 3))

cells_lipid_expr_1 %>%
  make_heatmap_jf_e1(design_expr_1, gr_expr_1, contrasts_expr_1, "Cells Lipids in Experiment 1", 0.01, 
                     seq(-2, 2, length.out = 100), -2:2, -2:2)

media_lipid_expr_1 %>% select(-blank_solvent) %>%
  make_heatmap_jf_e1(design_expr_1, gr_expr_1, contrasts_expr_1, "Media Lipids in Experiment 1", 0.01,
                     seq(-1, 1, length.out = 100), -1:1, -1:1)

cells_metab_expr_1 %>%
  make_heatmap_jf_e1(design_expr_1, gr_expr_1, contrasts_expr_1, "Cells Metabolites in Experiment 1", 0.01,
                     seq(-2, 2, length.out = 100), -2:2, -2:2)

media_metab_expr_1 %>% select(-solvent_blank_media) %>%
  make_heatmap_jf_e1(design_expr_1, gr_expr_1, contrasts_expr_1, "Media Metabolites in Experiment 1", 0.01,
                     seq(-1.3, 1.3, length.out = 100), seq(-1.3, 1.3, length.out = 3), seq(-1.3, 1.3, length.out = 3))

cells_lipid_expr_1 %>%
  make_heatmap_jf_e1(design_expr_1, gr_expr_1, contrasts_expr_1, "Cells Lipids in Experiment 1", 0.1, 
                     seq(-1, 1, length.out = 100), -1:1, -1:1)

media_lipid_expr_1 %>% select(-blank_solvent) %>%
  make_heatmap_jf_e1(design_expr_1, gr_expr_1, contrasts_expr_1, "Media Lipids in Experiment 1", 0.1,
                     seq(-1, 1, length.out = 100), -1:1, -1:1)

cells_metab_expr_1 %>%
  make_heatmap_jf_e1(design_expr_1, gr_expr_1, contrasts_expr_1, "Cells Metabolites in Experiment 1", 0.1,
                     seq(-2, 2, length.out = 100), -2:2, -2:2)

media_metab_expr_1 %>% select(-solvent_blank_media) %>%
  make_heatmap_jf_e1(design_expr_1, gr_expr_1, contrasts_expr_1, "Media Metabolites in Experiment 1", 0.1,
                     seq(-2, 2, length.out = 100), -2:2, -2:2)
