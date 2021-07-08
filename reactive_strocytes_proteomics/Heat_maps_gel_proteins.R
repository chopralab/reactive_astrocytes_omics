library(plyr)
library(tidyverse)
library(RColorBrewer)
library(RColorBrewer)

df <- read.csv('4x_full_for_analysis.csv')
df[, 8:19][df[, 8:19] == 0] <- 1
# df_sig_p <- filter(df, df$FDR<0.1)


#-average_C_lfq, -average_R_lfq,

df %>% mutate(proteins = make.unique(Protein.ID)) %>%
  filter(FDR<0.1) %>%
  select(-Majority.protein.IDs, -logFC, -PValue, -Fold_change, 
         -MS.MS.count.C1, -MS.MS.count.C2, -MS.MS.count.C3, -MS.MS.count.C4, -MS.MS.count.C5,
         -MS_C_average, -MS.MS.count.R1, -MS.MS.count.R2, -MS.MS.count.R3, -MS.MS.count.R4,
         -MS.MS.count.R5, -MS_R_average, -FDR, -Gene.names,
         -C1, -C2, -C3, -C4, -C5, -R1, -R2, -R3, -R4, -R5) %>%
  rowwise() %>%
  mutate(mean = mean(c(average_C_lfq, average_R_lfq))) %>%
  mutate_if(is.numeric, log10) %>%
  mutate_if(is.numeric, list(~ . - mean)) %>%
  mutate(sd = sd(c(average_C_lfq, average_R_lfq))) %>%
  mutate_if(is.numeric, list(~ . / sd)) %>%
  select(-Protein.ID, -mean, -sd) ->
  vals

vals %>%
  column_to_rownames("proteins") %>%
  as.matrix() %>%
  pheatmap::pheatmap(color = colorRampPalette(rev(brewer.pal(n = 7, name =
                                                               "RdBu")))(100), 
                     main = "Astrocyte Cells Proteins", 
                     cluster_rows = TRUE,
                     cluster_cols = FALSE, breaks = seq(-0.2, 0.2, length.out = 100),legend_breaks = -0.2:0.2, legend_labels = -0.2:0.2,
                     fontsize = 10, fontsize_row = 6)









# make_heatmap_jf_e1 <- function(tp, design_mat, gr, contrasts, title = "Heat-map", FDR = 0.10, breaks, list1, list2) {
# 
#   tp %>%
#     mutate(lipid = make.unique(lipid)) %>%
#     filter(lipid %in% DElist) %>%
#     select(-Transition, -type) %>%
#     rowwise() %>%
#     mutate(mean = mean(c(C1, C2, C3, C4, C5, R1, R2, R3, R4, R5))) %>%
#     mutate_if(is.numeric, log2) %>%
#     mutate_if(is.numeric, list(~ . - mean)) %>%
#     select(-blank, -mean) ->
#     vals
#   
#   vals %>%
#     column_to_rownames("lipid") %>%
#     as.matrix() %>%
#     pheatmap::pheatmap(color = colorRampPalette(rev(brewer.pal(n = 7, name =
#                                                                  "RdBu")))(100), 
#                        main = title, 
#                        cluster_rows = FALSE,
#                        cluster_cols = FALSE, breaks = breaks, legend_breaks = list1, legend_labels = list2,
#                        fontsize = 10, fontsize_row = 6)
# }
# 
# cells_lipid_expr_1 %>%
#   make_heatmap_jf_e1(design_expr_1, gr_expr_1, contrasts_expr_1, "Cells Lipids in Experiment 1", 0.001, 
#                      seq(-2, 2, length.out = 100), -2:2, -2:2)
