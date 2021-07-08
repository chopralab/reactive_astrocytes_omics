library(ggplot2)
library(ggrepel)
library(tidyverse)

df <- read.csv(file = '4x_full_for_analysis_avg_only.csv')
df[, 6:9][df[, 6:9] == 0] <- 1


df_sig_p <- filter(df, df$FDR<0.1)


mid_logfc <- filter(df, 2 > df$logFC & df$logFC > -2)

high_logfc <- filter(df, 2 < df$logFC)
high_logfc_labelling <- filter(df, 1 < df$logFC & df$logFC < 2)

low_logfc <- filter(df, df$logFC < -2)

low_and_high <- rbind(high_logfc, low_logfc)

# top97 <- filter(df, df$FDR<2.640000e-26)

### changing zeros in column 6 and 7 to 1's

mid<-5e+09

## Normal scale
df %>% 
  ggplot(aes(x=average_C_lfq,y=average_R_lfq, color=Fold_change)) + 
  geom_point(alpha=1, color="#595855", size = 3) +
  
  geom_point(data=df_sig_p, 
             aes(x=average_C_lfq,y=average_R_lfq),
             size=3) + scale_color_gradient2(midpoint=mid, low="red",mid="#DC9313",
                                             high="blue",  guide = "colourbar" )+
  # annotate("segment", x = 0, xend = 2e+10, y = 1285940000, yend = 0.7e+10, colour = "black")+
  # annotate("text",label="C3", x=2e+10, y=0.7e+10, size=5,fontface="bold.italic", color="black")+
  # 
  # annotate("segment", x = 0, xend = 5e+9, y = 171364000, yend = 2e+10, colour = "black")+
  # annotate("text",label="Serping1", x=5e+9, y=2e+10, size=5,fontface="bold.italic", color="black")+
  # 
  # annotate("segment", x = 615012000, xend = 4e+10, y = 353988000, yend =1e+10, colour = "black")+
  # annotate("text",label="Apoe", x=4e+10, y=1e+10, size=5,fontface="bold.italic", color="black")+  	 
  # 
  # annotate("segment", x = 236912000, xend = 3e+10, y = 163656000, yend = 1.5e+10, colour = "black")+
  # annotate("text",label="Sparc", x=3e+10, y=1.5e+10, size=5,fontface="bold.italic", color="black")+	 
  # 
  # annotate("segment", x = 118808000, xend = 1.5e+10, y = 39262600, yend = 3e+10, colour = "black")+
  # annotate("text",label="Cst3", x=1.5e+10, y=3e+10, size=5,fontface="bold.italic", color="black")+  	 
  # 
  #   
  # annotate("segment", x = 36163200, xend = 4e+10, y = 167992000, yend = 2.3e+9, colour = "black")+
  # annotate("text",label="Clu", x=4e+10, y=2.3e+9, size=5,fontface="bold.italic", color="black")+
  
  ggtitle("Astrocyte Gel Proteins") +
  theme(plot.title = element_text(hjust = 0.5))+
  ggrepel::geom_text_repel(data = top100, aes(label = Gene.names))

##Log10 scale
mid<-0

df %>% 
  ggplot(aes(x=log10(MS_C_average),y=log10(MS_R_average), color=logFC)) + 
  
  geom_point(data = low_logfc, alpha=1, color="blue", size = 1.5) +
  
  geom_point(data = high_logfc, alpha=1, color="red", size = 1.5) +
  
  geom_point(data=mid_logfc,aes(x=log10(MS_C_average),y=log10(MS_R_average)),size=1.5)+
  
  scale_color_gradient2(midpoint=mid, low="blue",mid="#DC9313",
                                             high="red", space ="Lab", limits = c(-2, 2), 
                                             breaks=c(-2, 0, 2.))+
  theme_bw()+
  # annotate("segment", x = log10(1), xend = 2.5, y = log10(1285940000), yend = 10, colour = "black")+
  # annotate("text",label="C3", x=2.5, y=10, size=5,fontface="bold.italic", color="black")+
  # 
  # annotate("segment", x = log10(1), xend = 2.5, y = log10(171364000), yend = 7, colour = "black")+
  # annotate("text",label="Serping1", x=2.5, y=7, size=5,fontface="bold.italic", color="black")+
  # 
  # annotate("segment", x = log10(615012000), xend = 11, y = log10(353988000), yend =8, colour = "black")+
  # annotate("text",label="Apoe", x=11, y=8, size=5,fontface="bold.italic", color="black")+
  # 
  # annotate("segment", x = log10(236912000), xend = 7, y = log10(163656000), yend = 11, colour = "black")+
  # annotate("text",label="Sparc", x=7, y=11, size=5,fontface="bold.italic", color="black")+
  # 
  # annotate("segment", x = log10(118808000), xend = 10, y = log10(39262600), yend = 6, colour = "black")+
  # annotate("text",label="Cst3", x=10, y=6, size=5,fontface="bold.italic", color="black")+
  # 
  # 
  # annotate("segment", x = log10(36163200), xend = 4.5, y = log10(167992000), yend = 9, colour = "black")+
  # annotate("text",label="Clu", x= 4.5, y=9, size=5,fontface="bold.italic", color="black")+
  
  ggtitle("Astrocyte Cell Proteins") +
  theme(plot.title = element_text(hjust = 0.5))+
  ggrepel::geom_text_repel(data = high_logfc_labelling, aes(label = Gene.names))
  
  

