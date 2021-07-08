library(ggplot2)
library(ggrepel)
library(tidyverse)

df <- read.csv(file = 'EKO_R_VS_WT_R_cells_full.csv')

df_sig_p <- filter(df, df$FDR<0.1)

##Log10 scale
mid<-0

df %>% 
  ggplot(aes(x=log10(EKO_average),y=log10(WT_average), color=logFC)) + 
  geom_point(alpha=1, color="grey", size = 1.5) +
  
  geom_point(data=df_sig_p, 
             aes(x=log10(EKO_average),y=log10(WT_average)),
             size=1.5) + scale_color_gradient2(midpoint=mid, low="blue",mid="#DC9313",
                                             high="red", space ="Lab", limits = c(-2.3, 2.3), 
                                             breaks=c(-2, 0, 2.))+
  theme_bw()+
  ggtitle("EKO-R Vs WT-R - Cells") +
  theme(plot.title = element_text(hjust = 0.5))+
  ggrepel::geom_text_repel(data = df_sig_p, aes(label = lipid))
  

########## MEDIA ###################

df_media <- read.csv(file = 'EKO_R_VS_WT_R_media_full.csv')
df_sig_p_media <- filter(df_media, df_media$FDR<0.1)

##Log10 scale
mid<-0

df_media %>% 
  ggplot(aes(x=log10(EKO_average),y=log10(WT_average), color=logFC)) + 
  geom_point(alpha=1, color="grey", size = 1.5) +
  
  geom_point(data=df_sig_p_media, 
             aes(x=log10(EKO_average),y=log10(WT_average)),
             size=1.5) + scale_color_gradient2(midpoint=mid, low="blue",mid="#DC9313",
                                               high="red", space ="Lab", limits = c(-1, 1), 
                                               breaks=c(-1, 0, 1.))+
  theme_bw()+
  ggtitle("EKO-R Vs WT-R - Media") +
  theme(plot.title = element_text(hjust = 0.5))+
  ggrepel::geom_text_repel(data = df_sig_p_media, aes(label = lipid))
