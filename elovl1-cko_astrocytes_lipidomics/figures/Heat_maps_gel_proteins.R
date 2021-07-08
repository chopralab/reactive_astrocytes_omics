library(plyr)
library(tidyverse)
library(RColorBrewer)

df <- read.csv('EKO_R_VS_WT_R_cells_full.csv')


#-average_C_lfq, -average_R_lfq,

df %>% mutate(lipids = make.unique(lipid)) %>%
  filter(FDR<0.1) %>%
  select( -logFC, -logCPM,	-LR, -PValue, -FDR,-Transition,	-type,
          -EKO_R1_Cells_r1,	-EKO_R1_Cells_r2,	-EKO_R1_Cells_r3,	-EKO_R2_Cells_r1,	
          -EKO_R2_Cells_r2,	-EKO_R2_Cells_r3,	-EKO_R3_Cells_r1,	-EKO_R3_Cells_r2,	
          -EKO_R3_Cells_r3,	-WT_R1_Cells_r1,	-WT_R1_Cells_r2,	-WT_R1_Cells_r3,
          -WT_R2_Cells_r1,	-WT_R2_Cells_r2,	-WT_R2_Cells_r3,	-WT_R3_Cells_r1,	
          -WT_R3_Cells_r2,	-WT_R3_Cells_r3,
          -blank) %>%
  rowwise() %>%
  mutate(mean = mean(c(EKO_average, WT_average))) %>%
  mutate_if(is.numeric, log10) %>%
  mutate_if(is.numeric, list(~ . - mean)) %>%
  mutate(sd = sd(c(EKO_average, WT_average))) %>%
  mutate_if(is.numeric, list(~ . / sd)) %>%
  select(-lipid, -mean, -sd) ->
  vals

vals %>%
  column_to_rownames("lipids") %>%
  as.matrix() %>%
  pheatmap::pheatmap(color = colorRampPalette(rev(brewer.pal(n = 7, name =
                                                               "RdBu")))(100), 
                     main = "EKO-R Vs WT-R Cells", 
                     cluster_rows = TRUE,
                     cluster_cols = FALSE, breaks = seq(-2.2, 1.5, length.out = 100),
                     legend_breaks = c(-2, -1, 0, 1), legend_labels = c(-2, -1, 0, 1),
                     fontsize = 10, fontsize_row = 6)

##############  MEDIA ################################

df <- read.csv('EKO_R_VS_WT_R_media_full.csv')


#-average_C_lfq, -average_R_lfq,

df %>% mutate(lipids = make.unique(lipid)) %>%
  filter(FDR<0.1) %>%
  select( -logFC, -logCPM,	-LR, -PValue, -FDR,-Transition,	-type,
          -EKO.R1.Media_r1, -EKO.R1.Media_r3, -EKO.R2.Media_r1, -EKO.R2.Media_r2, -EKO.R2.Media_r3,
          -EKO.R3.Media_r1, -EKO.R3.Media_r2, -EKO.R3.Media_r3, -WT.R1.Media_r1, -WT.R1.Media_r3,
          -WT.R2.Media_r1, -WT.R2.Media_r1, -WT.R2.Media_r2, -WT.R2.Media_r3, -WT.R3.Media_r1,
          -WT.R3.Media_r2, -WT.R3.Media_r3,
          -blank) %>%
  rowwise() %>%
  mutate(mean = mean(c(EKO_average, WT_average))) %>%
  mutate_if(is.numeric, log10) %>%
  mutate_if(is.numeric, list(~ . - mean)) %>%
  mutate(sd = sd(c(EKO_average, WT_average))) %>%
  mutate_if(is.numeric, list(~ . / sd)) %>%
  select(-lipid, -mean, -sd) ->
  vals

vals %>%
  column_to_rownames("lipids") %>%
  as.matrix() %>%
  pheatmap::pheatmap(color = colorRampPalette(rev(brewer.pal(n = 7, name =
                                                               "RdBu")))(100), 
                     main = "EKO-R Vs WT-R Media", 
                     cluster_rows = TRUE,
                     cluster_cols = FALSE, breaks = seq(-2, 1.5, length.out = 100),
                     legend_breaks = c(-2, -1, 0, 1), legend_labels = c(-2, -1, 0, 1),
                     fontsize = 10, fontsize_row = 6)



