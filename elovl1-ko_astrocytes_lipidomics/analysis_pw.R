library(plyr)
library(tidyverse)

library(readxl)

experiment_1_names <- c("c1", "c2", "c3",
                        "R1", "R2", "R3",
                        "blank")

# experiment_2_names <- c("Control1", "Control2", "Control3", "Control4",
#                         "RAC1", "RAC2", "RAC3", "RAC4",
#                         "blank")

lipid_types <- c("FFA" = "FFA", "PC" = "PC")

# metob_types <- c("NEG1" = "NEG1", "NEG2" = "NEG2", "POS1" = "POS1", "POS2" = "POS2")

read_lipids <- function(prefix, marker, exp_names, types, add_blank = NULL) {
  lapply(names(types), function(type){
    
    dir_name = paste0(prefix, type, "/")
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
  "raw_data1/cells_lipids1/", "cells",
  experiment_1_names, lipid_types
)

media_lipid_expr_1 = read_lipids(
  "raw_data1/media_lipids1/", "media", 
  experiment_1_names, lipid_types
)


library(edgeR)

# Groups for experiment 1
gr_expr_1 = c("C", "C", "C",
              "R", "R", "R",
              "Blank") %>%
  factor(levels = c("Blank", "C", "R"))

design_expr_1 = model.matrix(~gr_expr_1)

contrasts_expr_1 = makeContrasts(
  H = gr_expr_1R - gr_expr_1C,
  levels = design_expr_1
)

# # Groups for experiment 2
# gr_expr_2 = c("C", "C", "C", "C",
#               "R", "R", "R", "R",
#               "Blank") %>%
#   factor(levels = c("Blank", "C", "R"))
# 
# design_expr_2 = model.matrix(~gr_expr_2)
# 
# contrasts_expr_2 = makeContrasts(
#   H = gr_expr_2R - gr_expr_2C,
#   levels = design_expr_2
# )

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
  # select(-blank_solvent) %>%
  experiment_1_helper

# cm_e1_tbl <-
#   cells_metab_expr_1 %>%
#   experiment_1_helper
# 
# mm_e1_tbl <-
#   media_metab_expr_1 %>%
#   # select(-solvent_blank_media) %>%
#   experiment_1_helper

cl_e1_tbl

ml_e1_tbl

# cm_e1_tbl
# 
# mm_e1_tbl

# experiment_2_helper <- function(df) {
#   df %>%
#     perform_analysis_raw(design_expr_2, gr_expr_2) %>%
#     calculate_significance(contrasts_expr_2) %>%
#     topTags(1500000) %>%
#     as.data.frame() %>%
#     rownames_to_column("lipid") %>%
#     as_tibble()
#   
# }
# 
# cl_e2_tbl <-
#   cells_lipid_expr_2 %>%
#   experiment_2_helper
# 
# ml_e2_tbl <-
#   media_lipid_expr_2 %>%
#   select(-blank_solvent) %>%
#   experiment_2_helper
# 
# cm_e2_tbl <-
#   cells_metab_expr_2 %>%
#   experiment_2_helper
# 
# mm_e2_tbl <-
#   media_metab_expr_2 %>%
#   select(-solvent_blank_media) %>%
#   experiment_2_helper
# 
# cl_e2_tbl
# 
# ml_e2_tbl
# 
# cm_e2_tbl
# 
# mm_e2_tbl

make_volcano_plot <- function(df, title) {
  df %>%
    mutate(sig = factor(FDR < 0.10)) %>%
    ggplot(aes(logFC, -log10(FDR), color = sig)) +
    geom_point() +
    scale_color_manual(values = c("FALSE" = "black", "TRUE" = "red")) +
    guides(color = F) +
    ggtitle(title)
}

cl_e1_tbl %>% make_volcano_plot("Cells Lipids Experiment 1")

ml_e1_tbl %>% make_volcano_plot("Media Lipids Experiment 1")

# cm_e1_tbl %>% make_volcano_plot("Cells Metabolites Experiment 1")
# 
# mm_e1_tbl %>% make_volcano_plot("Media Metabolites Experiments 1")

# cl_e2_tbl %>% make_volcano_plot("Cells Lipids Experiment 2")
# 
# ml_e2_tbl %>% make_volcano_plot("Media Lipids Experiment 2")
# 
# cm_e2_tbl %>% make_volcano_plot("Cells Metabolites Experiment 2")
# 
# mm_e2_tbl %>% make_volcano_plot("Media Metabolites Experiments 2")


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

cl_e1_tbl %>% write_summary_and_results(cells_lipid_expr_1, "cl_e1")
ml_e1_tbl %>% write_summary_and_results(media_lipid_expr_1, "ml_e1")
# cm_e1_tbl %>% write_summary_and_results(cells_metab_expr_1, "cm_e1")
# mm_e1_tbl %>% write_summary_and_results(media_metab_expr_1, "mm_e1")

# cl_e2_tbl %>% write_summary_and_results(cells_lipid_expr_2, "cl_e2")
# ml_e2_tbl %>% write_summary_and_results(media_lipid_expr_2, "ml_e2")
# cm_e2_tbl %>% write_summary_and_results(cells_metab_expr_2, "cm_e2")
# mm_e2_tbl %>% write_summary_and_results(media_metab_expr_2, "mm_e2")


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
    ) + xlim(-2,2) + ylim(-2.5,2.5)+
    ggtitle(title) +
    cowplot::theme_cowplot()
}

cells_lipid_expr_1 %>%
  make_pca_plot(design_expr_1, gr_expr_1, contrasts_expr_1, "Cells Lipids in Experiment 1")

media_lipid_expr_1 %>% #select(-blank_solvent) %>%
  make_pca_plot(design_expr_1, gr_expr_1, contrasts_expr_1, "Media Lipids in Experiment 1")
# 
# cells_metab_expr_1 %>%
#   make_pca_plot(design_expr_1, gr_expr_1, contrasts_expr_1, "Cells Metabolites in Experiment 1")
# 
# media_metab_expr_1 %>% select(-solvent_blank_media) %>%
#   make_pca_plot(design_expr_1, gr_expr_1, contrasts_expr_1, "Media Metabolites in Experiment 1")

# cells_lipid_expr_2 %>%
#   make_pca_plot(design_expr_2, gr_expr_2, contrasts_expr_2, "Cells Lipids in Experiment 2")
# 
# media_lipid_expr_2 %>% select(-blank_solvent) %>%
#   make_pca_plot(design_expr_2, gr_expr_2, contrasts_expr_2, "Media Lipids in Experiment 2")
# 
# cells_metab_expr_2 %>%
#   make_pca_plot(design_expr_2, gr_expr_2, contrasts_expr_2, "Cells Metabolites in Experiment 2")
# 
# media_metab_expr_2 %>% select(-solvent_blank_media) %>%
#   make_pca_plot(design_expr_2, gr_expr_2, contrasts_expr_2, "Media Metabolites in Experiment 2")

#heatmaps

make_heatmap <- function(tp, design_mat, gr, contrasts, title = "Heat-map") {
  
  DElist <-
    tp %>%
    get_DE_lipids(design_mat, gr, contrasts)
  
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
    mutate_if(is.numeric, log2) %>%
    mutate_if(is.numeric, list(~ . - blank)) %>%
    select(-blank) %>%
    column_to_rownames("lipid") %>%
    as.matrix() %>%
    pheatmap::pheatmap(main = title)
}

cells_lipid_expr_1 %>%
  make_heatmap(design_expr_1, gr_expr_1, contrasts_expr_1, "Cells Lipids in Experiment 1")

media_lipid_expr_1 %>% #select(-blank_solvent) %>%
  make_heatmap(design_expr_1, gr_expr_1, contrasts_expr_1, "Media Lipids in Experiment 1")

# cells_metab_expr_1 %>%
#   make_heatmap(design_expr_1, gr_expr_1, contrasts_expr_1, "Cells Metabolites in Experiment 1")
# 
# 
# media_metab_expr_1 %>% select(-solvent_blank_media) %>%
#   make_heatmap(design_expr_1, gr_expr_1, contrasts_expr_1, "Media Metabolites in Experiment 1")


# cells_lipid_expr_2 %>%
#   make_heatmap(design_expr_2, gr_expr_2, contrasts_expr_2, "Cells Lipids in Experiment 2")
# 
# media_lipid_expr_2 %>% select(-blank_solvent) %>%
#   make_heatmap(design_expr_2, gr_expr_2, contrasts_expr_2, "Media Lipids in Experiment 2")
# 
# cells_metab_expr_2 %>%
#   make_heatmap(design_expr_2, gr_expr_2, contrasts_expr_2, "Cells Metabolites in Experiment 2")
# 
# media_metab_expr_2 %>% select(-solvent_blank_media) %>%
#   make_heatmap(design_expr_2, gr_expr_2, contrasts_expr_2, "Media Metabolites in Experiment 2")



make_heatmap_jf_e1 <- function(tp, design_mat, gr, contrasts, title = "Heat-map", FDR = 0.10) {
  
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
    pheatmap::pheatmap(main = title)
}

cells_lipid_expr_1 %>%
  make_heatmap_jf_e1(design_expr_1, gr_expr_1, contrasts_expr_1, "Cells Lipids in Experiment 1", 0.001)

media_lipid_expr_1 %>% select(-blank_solvent) %>%
  make_heatmap_jf_e1(design_expr_1, gr_expr_1, contrasts_expr_1, "Media Lipids in Experiment 1", 0.001)

# cells_metab_expr_1 %>%
#   make_heatmap_jf_e1(design_expr_1, gr_expr_1, contrasts_expr_1, "Cells Metabolites in Experiment 1", 0.001)
# 
# media_metab_expr_1 %>% select(-solvent_blank_media) %>%
#   make_heatmap_jf_e1(design_expr_1, gr_expr_1, contrasts_expr_1, "Media Metabolites in Experiment 1", 0.001)

# make_heatmap_jf_e2 <- function(tp, design_mat, gr, contrasts, title = "Heat-map", FDR = 0.10) {
#   
#   DElist <-
#     tp %>%
#     get_DE_lipids(design_mat, gr, contrasts, FDR)
#   
#   if(length(DElist) == 0) {
#     cat("No significant lipids for ", title)
#     return()
#   }
#   
#   if(length(DElist) == 1) {
#     cat("Single significant lipids for ", title, " is ", DElist[1])
#     return()
#   }
#   
#   tp %>%
#     mutate(lipid = make.unique(lipid)) %>%
#     filter(lipid %in% DElist) %>%
#     select(-Transition, -type) %>%
#     rowwise() %>%
#     mutate(mean = mean(c(Control1, Control2, Control3, Control4,
#                          RAC1, RAC2, RAC3, RAC4))) %>%
#     mutate_if(is.numeric, log2) %>%
#     mutate_if(is.numeric, list(~ . - mean)) %>%
#     select(-blank, -mean) ->
#     vals
#   
#   vals %>%
#     column_to_rownames("lipid") %>%
#     as.matrix() %>%
#     pheatmap::pheatmap(main = title)
# }
# 
# cells_lipid_expr_2 %>%
#   make_heatmap_jf_e2(design_expr_2, gr_expr_2, contrasts_expr_2, "Cells Lipids in Experiment 2", 0.001)
# 
# media_lipid_expr_2 %>% select(-blank_solvent) %>%
#   make_heatmap_jf_e2(design_expr_2, gr_expr_2, contrasts_expr_2, "Media Lipids in Experiment 2", 0.001)
# 
# cells_metab_expr_2 %>%
#   make_heatmap_jf_e2(design_expr_2, gr_expr_2, contrasts_expr_2, "Cells Metabolites in Experiment 2", 0.001)
# 
# media_metab_expr_2 %>% select(-solvent_blank_media) %>%
#   make_heatmap_jf_e2(design_expr_2, gr_expr_2, contrasts_expr_2, "Media Metabolites in Experiment 2", 0.001)

write_summary_and_results_combined_pp <- function(tbl1, tbl2, df1, df2, name, p.value = 0.1) {
  
  tbl1 %>%
    merge(df1) %>%
    as_tibble() %>%
    filter(FDR < p.value) -> results1
  
  tbl2 %>%
    merge(df2) %>%
    as_tibble() %>%
    filter(FDR < p.value) -> results2
  
  merge(results1, results2, by="lipid", suffixes = c(".cells", ".media")) %>%
    as_tibble() -> cl.ml
  
  cl.ml %>% arrange(FDR.cells, FDR.media)
}

write_summary_and_results_combined_pp(cl_e1_tbl, ml_e1_tbl,
                                      cells_lipid_expr_1, media_lipid_expr_1) #%>%

#write_csv("results/cl_w_ml_e1.csv")