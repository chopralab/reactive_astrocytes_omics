## For Cells Lipids



library(ggplot2)
library(ggrepel)
library(tidyverse)

# with pvalue < 0.05

df <- read.csv(file = 'All_absolute_intensities_foldchange_averages.csv')

df_sig_p <- filter(df, df$PValue<0.05)

mid<-mean(df_sig_p$Fold.Change)

df %>% 
  ggplot(aes(x=averageC,y=averageR, color=Fold.Change)) + 
  geom_point(alpha=1, color="darkgrey", size = 3) +

  geom_point(data=df_sig_p, 
             aes(x=averageC,y=averageR),
             size=3) + scale_color_gradient2(midpoint=mid, low="red", mid="white",
                                             high="blue", space ="Lab" )+
  ggrepel::geom_text_repel(data = df_sig_p, aes(label = df_sig_p$ï..Lipid.name))
  # annotate("text",label="16:0 Cholesteryl ester", x=4405.840315, y=4292.752293, size=5,fontface="bold.italic", color="black")



########################### Based on FDR < 0.1 ###########################################

library(ggplot2)
library(ggrepel)

df <- read.csv(file = 'All_absolute_intensities_foldchange_averages.csv')

df_sig_p <- filter(df, df$FDR<0.1)
PC_FFA <- filter(df_sig_p, df_sig_p$lipid_class=='PC_SM' | df_sig_p$lipid_class=='FFA') ## Just for naming

mid<-mean(df_sig_p$logFC)

df %>% 
  ggplot(aes(x=averageC,y=averageR, color=logFC)) + 
  geom_point(alpha=1, color="#595855", size = 3) +
  
  geom_point(data=df_sig_p, 
             aes(x=averageC,y=averageR),
             size=3) + scale_color_gradient2(midpoint=mid, low="red",mid="#DC9313",
                                             high="blue", space ="Lab" )+
  ggrepel::geom_text_repel(data = PC_FFA, aes(label = ï..Lipid.name))
# annotate("text",label="16:0 Cholesteryl ester", x=4405.840315, y=4292.752293, size=5,fontface="bold.italic", color="black")

##Log10 scale
mid<- 0.0
# high<- max(df_sig_p$logFC)
# low<- min(df_sig_p$logFC)


df %>% 
  ggplot(aes(x=log10(averageC),y=log10(averageR), color=logFC)) + 
  geom_point(alpha=1, color="grey", size = 3) +
  
  geom_point(data=df_sig_p, 
             aes(x=log10(averageC),y=log10(averageR)),
             size=3) + scale_color_gradient2(midpoint=mid, low="Blue",mid="#DC9313",
                                             high="Red", space ="Lab", limits = c(-1., 1.), 
                                             breaks= c(-1.,-0.5,0, 0.5, 1.) )+theme_bw()+
  ggtitle("Astrocyte Cell Lipids") +
  theme(plot.title = element_text(hjust = 0.5)) +
  
  ggrepel::geom_text_repel(data = PC_FFA, aes(label = ï..Lipid.name))


##Only the ones with less than 1e+05 averageC and 125000 averageR

library(ggplot2)
library(ggrepel)

df <- read.csv(file = 'All_absolute_intensities_foldchange_averages.csv')

sig_groups_pvalue <- as.factor(df$sig_groups_pvalue)
sig_groups_fdr <- as.factor(df$sig_groups_fdr)

df_smallc <- filter(df, df$averageC<1e+05)
df_smallr <- filter(df_smallc, df_smallc$averageR<125000)

df_sig_p <- filter(df_smallr, df_smallr$FDR<0.1)

mid<-mean(df_sig_p$Fold.Change)

df_smallr %>% 
  ggplot(aes(x=averageC,y=averageR, color=Fold.Change)) + 
  geom_point(alpha=1, color="#595855", size = 3) +
  
  geom_point(data=df_sig_p, 
             aes(x=averageC,y=averageR),
             size=3) + scale_color_gradient2(midpoint=mid, low="red",mid="#DC9313",
                                             high="blue", space ="Lab" )



## Media Lipids

## Based on FDR < 0.1

library(ggplot2)
library(ggrepel)

df <- read.csv(file = 'media_lipids_data.csv')


df_sig_p <- filter(df, df$FDR<0.1)
PC_FFA <- filter(df_sig_p, df_sig_p$lipid_class=='PC_SM' | df_sig_p$lipid_class=='FFA') ## Just for naming

b <- c(-1.5, 0, 1)
mid<-0

df %>% 
  ggplot(aes(x=AverageC,y=AverageR, color=logFC)) + 
  geom_point(alpha=1, color="grey", size = 3) +
  
  geom_point(data=df_sig_p, 
             aes(x=AverageC,y=AverageR),
             size=3) +
  scale_color_gradient2(midpoint=mid, low="blue",mid="#DC9313",
                                             high="red", space ="Lab", limits = c(-1.5,1.5))+
  
  ggrepel::geom_text_repel(data = PC_FFA, aes(label = ï..Lipid.name))
# annotate("text",label="16:0 Cholesteryl ester", x=4405.840315, y=4292.752293, size=5,fontface="bold.italic", color="black")

### LOG10 scale

df %>% 
  ggplot(aes(x=log10(AverageC),y=log10(AverageR), color=logFC)) + 
  geom_point(alpha=1, color="grey", size = 3) +
  
  geom_point(data=df_sig_p, 
             aes(x=log10(AverageC),y=log10(AverageR)),
             size=3) + scale_color_gradient2(midpoint=mid, low="blue",mid="#DC9313",
                                             high="red", space ="Lab", limits = c(-1,1), 
                                             breaks= c(-1.,-0.5,0, 0.5, 1.))+theme_bw()+
  ggtitle("Astrocyte Media Lipids") +
  theme(plot.title = element_text(hjust = 0.5)) + 
  ggrepel::geom_text_repel(data = PC_FFA, aes(label = ï..Lipid.name))



##Only the ones with 1e+05 averageC and 125000 averageR

library(ggplot2)
library(ggrepel)

df <- read.csv(file = 'media_lipids_data.csv')

df_smallc <- filter(df, df$AverageC<125000)
# df_smallr <- filter(df_smallc, df_smallc$AverageR<230000)

df_sig_p <- filter(df_smallc, df_smallc$FDR<0.1)

mid<-mean(df_sig_p$FoldChange)

df_smallc %>% 
  ggplot(aes(x=AverageC,y=AverageR, color=FoldChange)) + 
  geom_point(alpha=1, color="#595855", size = 3) +
  
  geom_point(data=df_sig_p, 
             aes(x=AverageC,y=AverageR),
             size=3) + scale_color_gradient2(midpoint=mid, low="red",mid="#DC9313",
                                             high="blue", space ="Lab" )


### plotting -log(Pvalue) vs. log2(foldchange) for media lipids



