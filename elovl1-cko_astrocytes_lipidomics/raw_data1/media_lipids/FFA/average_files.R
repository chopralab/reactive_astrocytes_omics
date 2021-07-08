library(plyr)
library(tidyverse)
library(readxl)
library(dplyr)

c1_r1_ffa <- read_excel("EKO-C1-Media_r1_FFA_results.xlsx")
c1_r2_ffa <- read_excel("EKO-C1-Media_r2_FFA_results.xlsx")
c1_r3_ffa <- read_excel("EKO-C1-Media_r3_FFA_results.xlsx")

c2_r1_ffa <- read_excel("EKO-C2-Media_r1_FFA_results.xlsx")
c2_r2_ffa <- read_excel("EKO-C2-Media_r2_FFA_results.xlsx")
c2_r3_ffa <- read_excel("EKO-C2-Media_r3_FFA_results.xlsx")


c3_r1_ffa <- read_excel("EKO-C3-Media_r1_FFA_results.xlsx")
c3_r2_ffa <- read_excel("EKO-C3-Media_r2_FFA_results.xlsx")
c3_r3_ffa <- read_excel("EKO-C3-Media_r3_FFA_results.xlsx")


R1_r1_ffa <- read_excel("EKO-R1-Media_r1_FFA_results.xlsx")
R1_r2_ffa <- read_excel("EKO-R1-Media_r2_FFA_results.xlsx")
R1_r3_ffa <- read_excel("EKO-R1-Media_r3_FFA_results.xlsx")

R2_r1_ffa <- read_excel("EKO-R2-Media_r1_FFA_results.xlsx")
R2_r2_ffa <- read_excel("EKO-R2-Media_r2_FFA_results.xlsx")
R2_r3_ffa <- read_excel("EKO-R2-Media_r3_FFA_results.xlsx")


R3_r1_ffa <- read_excel("EKO-R3-Media_r1_FFA_results.xlsx")
R3_r2_ffa <- read_excel("EKO-R3-Media_r2_FFA_results.xlsx")
R3_r3_ffa <- read_excel("EKO-R3-Media_r3_FFA_results.xlsx")

c1 <- as.data.frame(c(c1_r1_ffa[,3],c1_r2_ffa[,3],c1_r3_ffa[,3]), drop=FALSE)
c1<- data.frame(Means=rowMeans(c1[,-1]))
c1_1<- data.frame(c1_r1_ffa[,1:2])
C1_average <- bind_cols(c1_1,c1)

c2 <- as.data.frame(c(c2_r1_ffa[,3],c2_r2_ffa[,3],c2_r3_ffa[,3]), drop=FALSE)
c2<- data.frame(Means=rowMeans(c2[,-1]))
c2_1<- data.frame(c2_r1_ffa[,1:2])
C2_average <- bind_cols(c2_1,c2)

c3 <- as.data.frame(c(c3_r1_ffa[,3],c3_r2_ffa[,3],c3_r3_ffa[,3]), drop=FALSE)
c3<- data.frame(Means=rowMeans(c3[,-1]))
c3_1<- data.frame(c3_r1_ffa[,1:2])
C3_average <- bind_cols(c3_1,c3)

R1 <- as.data.frame(c(R1_r1_ffa[,3],R1_r2_ffa[,3],R1_r3_ffa[,3]), drop=FALSE)
R1<- data.frame(Means=rowMeans(R1[,-1]))
R1_1<- data.frame(R1_r1_ffa[,1:2])
R1_average <- bind_cols(R1_1,R1)

R2 <- as.data.frame(c(R2_r1_ffa[,3],R2_r3_ffa[,3]), drop=FALSE)
R2<- data.frame(Means=rowMeans(R2[,-1]))
R2_1<- data.frame(R2_r1_ffa[,1:2])
R2_average <- bind_cols(R2_1,R2)

R3 <- as.data.frame(c(R3_r1_ffa[,3],R3_r2_ffa[,3],R3_r3_ffa[,3]), drop=FALSE)
R3<- data.frame(Means=rowMeans(R3[,-1]))
R3_1<- data.frame(R3_r1_ffa[,1:2])
R3_average <- bind_cols(R3_1,R3)


library("writexl")

write_xlsx(C1_average, "C1_cells_FFA_results_E.xlsx")
write_xlsx(C2_average, "C2_cells_FFA_results_E.xlsx")
write_xlsx(C3_average, "C3_cells_FFA_results_E.xlsx")
write_xlsx(R1_average, "R1_cells_FFA_results_E.xlsx")
write_xlsx(R2_average, "R2_cells_FFA_results_E.xlsx")
write_xlsx(R3_average, "R3_cells_FFA_results_E.xlsx")

#################################################################################

c1_r1_ffa <- read_excel("WT-C1-Media_r1_FFA_results.xlsx")
c1_r2_ffa <- read_excel("WT-C1-Media_r2_FFA_results.xlsx")
c1_r3_ffa <- read_excel("WT-C1-Media_r3_FFA_results.xlsx")

c2_r1_ffa <- read_excel("WT-C2-Media_r1_FFA_results.xlsx")
c2_r2_ffa <- read_excel("WT-C2-Media_r2_FFA_results.xlsx")
c2_r3_ffa <- read_excel("WT-C2-Media_r3_FFA_results.xlsx")


c3_r1_ffa <- read_excel("WT-C3-Media_r1_FFA_results.xlsx")
c3_r2_ffa <- read_excel("WT-C3-Media_r2_FFA_results.xlsx")
c3_r3_ffa <- read_excel("WT-C3-Media_r3_FFA_results.xlsx")


R1_r1_ffa <- read_excel("WT-R1-Media_r1_FFA_results.xlsx")
R1_r2_ffa <- read_excel("WT-R1-Media_r2_FFA_results.xlsx")
R1_r3_ffa <- read_excel("WT-R1-Media_r3_FFA_results.xlsx")

R2_r1_ffa <- read_excel("WT-R2-Media_r1_FFA_results.xlsx")
R2_r2_ffa <- read_excel("WT-R2-Media_r2_FFA_results.xlsx")
R2_r3_ffa <- read_excel("WT-R2-Media_r3_FFA_results.xlsx")


R3_r1_ffa <- read_excel("WT-R3-Media_r1_FFA_results.xlsx")
R3_r2_ffa <- read_excel("WT-R3-Media_r2_FFA_results.xlsx")
R3_r3_ffa <- read_excel("WT-R3-Media_r3_FFA_results.xlsx")

c1 <- as.data.frame(c(c1_r1_ffa[,3],c1_r2_ffa[,3],c1_r3_ffa[,3]), drop=FALSE)
c1<- data.frame(Means=rowMeans(c1[,-1]))
c1_1<- data.frame(c1_r1_ffa[,1:2])
C1_average <- bind_cols(c1_1,c1)

c2 <- as.data.frame(c(c2_r1_ffa[,3],c2_r2_ffa[,3],c2_r3_ffa[,3]), drop=FALSE)
c2<- data.frame(Means=rowMeans(c2[,-1]))
c2_1<- data.frame(c2_r1_ffa[,1:2])
C2_average <- bind_cols(c2_1,c2)

c3 <- as.data.frame(c(c3_r1_ffa[,3],c3_r2_ffa[,3],c3_r3_ffa[,3]), drop=FALSE)
c3<- data.frame(Means=rowMeans(c3[,-1]))
c3_1<- data.frame(c3_r1_ffa[,1:2])
C3_average <- bind_cols(c3_1,c3)

R1 <- as.data.frame(c(R1_r1_ffa[,3],R1_r2_ffa[,3],R1_r3_ffa[,3]), drop=FALSE)
R1<- data.frame(Means=rowMeans(R1[,-1]))
R1_1<- data.frame(R1_r1_ffa[,1:2])
R1_average <- bind_cols(R1_1,R1)

R2 <- as.data.frame(c(R2_r1_ffa[,3],R2_r2_ffa[,3],R2_r3_ffa[,3]), drop=FALSE)
R2<- data.frame(Means=rowMeans(R2[,-1]))
R2_1<- data.frame(R2_r1_ffa[,1:2])
R2_average <- bind_cols(R2_1,R2)

R3 <- as.data.frame(c(R3_r1_ffa[,3],R3_r2_ffa[,3],R3_r3_ffa[,3]), drop=FALSE)
R3<- data.frame(Means=rowMeans(R3[,-1]))
R3_1<- data.frame(R3_r1_ffa[,1:2])
R3_average <- bind_cols(R3_1,R3)


library("writexl")

write_xlsx(C1_average, "C1_cells_FFA_results_WT.xlsx")
write_xlsx(C2_average, "C2_cells_FFA_results_WT.xlsx")
write_xlsx(C3_average, "C3_cells_FFA_results_WT.xlsx")
write_xlsx(R1_average, "R1_cells_FFA_results_WT.xlsx")
write_xlsx(R2_average, "R2_cells_FFA_results_WT.xlsx")
write_xlsx(R3_average, "R3_cells_FFA_results_WT.xlsx")










