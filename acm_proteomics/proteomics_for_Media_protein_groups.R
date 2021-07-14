# 
#  ggbiplot.r
#  
#  Copyright 2011 Vincent Q. Vu.
# 
#  This program is free software; you can redistribute it and/or
#  modify it under the terms of the GNU General Public License
#  as published by the Free Software Foundation; either version 2
#  of the License, or (at your option) any later version.
#  
#  This program is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#  
#  You should have received a copy of the GNU General Public License
#  along with this program; if not, write to the Free Software
#  Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301, USA.
#

#' Biplot for Principal Components using ggplot2
#'
#' @param pcobj           an object returned by prcomp() or princomp()
#' @param choices         which PCs to plot
#' @param scale           covariance biplot (scale = 1), form biplot (scale = 0). When scale = 1, the inner product between the variables approximates the covariance and the distance between the points approximates the Mahalanobis distance.
#' @param obs.scale       scale factor to apply to observations
#' @param var.scale       scale factor to apply to variables
#' @param pc.biplot       for compatibility with biplot.princomp()
#' @param groups          optional factor variable indicating the groups that the observations belong to. If provided the points will be colored according to groups
#' @param ellipse         draw a normal data ellipse for each group?
#' @param ellipse.prob    size of the ellipse in Normal probability
#' @param labels          optional vector of labels for the observations
#' @param labels.size     size of the text used for the labels
#' @param alpha           alpha transparency value for the points (0 = transparent, 1 = opaque)
#' @param circle          draw a correlation circle? (only applies when prcomp was called with scale = TRUE and when var.scale = 1)
#' @param var.axes        draw arrows for the variables?
#' @param varname.size    size of the text for variable names
#' @param varname.adjust  adjustment factor the placement of the variable names, >= 1 means farther from the arrow
#' @param varname.abbrev  whether or not to abbreviate the variable names
#'
#' @return                a ggplot2 plot
#' @export
#' @examples
#'   data(wine)
#'   wine.pca <- prcomp(wine, scale. = TRUE)
#'   print(ggbiplot(wine.pca, obs.scale = 1, var.scale = 1, groups = wine.class, ellipse = TRUE, circle = TRUE))
#'
ggbiplot <- function(pcobj, choices = 1:2, scale = 1, pc.biplot = TRUE, 
                     obs.scale = 1 - scale, var.scale = scale, 
                     groups = NULL, ellipse = FALSE, ellipse.prob = 0.68, 
                     labels = NULL, labels.size = 3, alpha = 1, 
                     var.axes = TRUE, 
                     circle = FALSE, circle.prob = 0.69, 
                     varname.size = 3, varname.adjust = 1.5, 
                     varname.abbrev = FALSE, ...)
{
  library(ggplot2)
  library(plyr)
  library(scales)
  library(grid)
  
  stopifnot(length(choices) == 2)
  
  # Recover the SVD
  if(inherits(pcobj, 'prcomp')){
    nobs.factor <- sqrt(nrow(pcobj$x) - 1)
    d <- pcobj$sdev
    u <- sweep(pcobj$x, 2, 1 / (d * nobs.factor), FUN = '*')
    v <- pcobj$rotation
  } else if(inherits(pcobj, 'princomp')) {
    nobs.factor <- sqrt(pcobj$n.obs)
    d <- pcobj$sdev
    u <- sweep(pcobj$scores, 2, 1 / (d * nobs.factor), FUN = '*')
    v <- pcobj$loadings
  } else if(inherits(pcobj, 'PCA')) {
    nobs.factor <- sqrt(nrow(pcobj$call$X))
    d <- unlist(sqrt(pcobj$eig)[1])
    u <- sweep(pcobj$ind$coord, 2, 1 / (d * nobs.factor), FUN = '*')
    v <- sweep(pcobj$var$coord,2,sqrt(pcobj$eig[1:ncol(pcobj$var$coord),1]),FUN="/")
  } else if(inherits(pcobj, "lda")) {
    nobs.factor <- sqrt(pcobj$N)
    d <- pcobj$svd
    u <- predict(pcobj)$x/nobs.factor
    v <- pcobj$scaling
    d.total <- sum(d^2)
  } else {
    stop('Expected a object of class prcomp, princomp, PCA, or lda')
  }
  
  # Scores
  choices <- pmin(choices, ncol(u))
  df.u <- as.data.frame(sweep(u[,choices], 2, d[choices]^obs.scale, FUN='*'))
  
  # Directions
  v <- sweep(v, 2, d^var.scale, FUN='*')
  df.v <- as.data.frame(v[, choices])
  
  names(df.u) <- c('xvar', 'yvar')
  names(df.v) <- names(df.u)
  
  if(pc.biplot) {
    df.u <- df.u * nobs.factor
  }
  
  # Scale the radius of the correlation circle so that it corresponds to 
  # a data ellipse for the standardized PC scores
  r <- sqrt(qchisq(circle.prob, df = 2)) * prod(colMeans(df.u^2))^(1/4)
  
  # Scale directions
  v.scale <- rowSums(v^2)
  df.v <- r * df.v / sqrt(max(v.scale))
  
  # Change the labels for the axes
  if(obs.scale == 0) {
    u.axis.labs <- paste('standardized PC', choices, sep='')
  } else {
    u.axis.labs <- paste('PC', choices, sep='')
  }
  
  # Append the proportion of explained variance to the axis labels
  u.axis.labs <- paste(u.axis.labs, 
                       sprintf('(%0.1f%% explained var.)', 
                               100 * pcobj$sdev[choices]^2/sum(pcobj$sdev^2)))
  
  # Score Labels
  if(!is.null(labels)) {
    df.u$labels <- labels
  }
  
  # Grouping variable
  if(!is.null(groups)) {
    df.u$groups <- groups
  }
  
  # Variable Names
  if(varname.abbrev) {
    df.v$varname <- abbreviate(rownames(v))
  } else {
    df.v$varname <- rownames(v)
  }
  
  # Variables for text label placement
  df.v$angle <- with(df.v, (180/pi) * atan(yvar / xvar))
  df.v$hjust = with(df.v, (1 - varname.adjust * sign(xvar)) / 2)
  
  # Base plot
  g <- ggplot(data = df.u, aes(x = xvar, y = yvar)) + 
    xlab(u.axis.labs[1]) + ylab(u.axis.labs[2]) + coord_equal()
  
  if(var.axes) {
    # Draw circle
    if(circle) 
    {
      theta <- c(seq(-pi, pi, length = 50), seq(pi, -pi, length = 50))
      circle <- data.frame(xvar = r * cos(theta), yvar = r * sin(theta))
      g <- g + geom_path(data = circle, color = muted('white'), 
                         size = 1/2, alpha = 1/3)
    }
    
    # Draw directions
    # g <- g +
    #   geom_segment(data = df.v,
    #                aes(x = 0, y = 0, xend = xvar, yend = yvar),
    #                arrow = arrow(length = unit(1/2, 'picas')), 
    #                color = muted('red'))
  }
  
  # Draw either labels or points
  if(!is.null(df.u$labels)) {
    if(!is.null(df.u$groups)) {
      g <- g + geom_text(aes(label = labels, color = groups), 
                         size = labels.size)
    } else {
      g <- g + geom_text(aes(label = labels), size = labels.size)      
    }
  } else {
    if(!is.null(df.u$groups)) {
      g <- g + geom_point(aes(color = groups), alpha = alpha)
    } else {
      g <- g + geom_point(alpha = alpha)      
    }
  }
  
  # Overlay a concentration ellipse if there are groups
  if(!is.null(df.u$groups) && ellipse) {
    theta <- c(seq(-pi, pi, length = 50), seq(pi, -pi, length = 50))
    circle <- cbind(cos(theta), sin(theta))
    
    ell <- ddply(df.u, 'groups', function(x) {
      if(nrow(x) <= 2) {
        return(NULL)
      }
      sigma <- var(cbind(x$xvar, x$yvar))
      mu <- c(mean(x$xvar), mean(x$yvar))
      ed <- sqrt(qchisq(ellipse.prob, df = 2))
      data.frame(sweep(circle %*% chol(sigma) * ed, 2, mu, FUN = '+'), 
                 groups = x$groups[1])
    })
    names(ell)[1:2] <- c('xvar', 'yvar')
    g <- g + geom_path(data = ell, aes(color = groups, group = groups))
  }
  
  # Label the variable axes
  # if(var.axes) {
  #   g <- g + 
  #     geom_text(data = df.v, 
  #               aes(label = varname, x = xvar, y = yvar, 
  #                   angle = angle, hjust = hjust), 
  #               color = 'darkred', size = varname.size)
  # }
  
  # if(var.axes) {
  #   g <- g +
  #     geom_point(data = df.v,
  #                color = 'darkred', size = 1)
  # }
  
  
  # Change the name of the legend for groups
  # if(!is.null(groups)) {
  #   g <- g + scale_color_brewer(name = deparse(substitute(groups)), 
  #                               palette = 'Dark2')
  # }
  
  # TODO: Add a second set of axes
  
  return(g)
}


##############################################################################

  
  # BV2 Proteomics
  
  ## Load data
  
library(plyr)
library(tidyverse)
protein_groups <- read_csv("Media_proteinGroups.csv")
protein_groups_lfq <- protein_groups %>%
  select(protein_ids = 'Majority protein IDs',
         R1 = 'LFQ intensity R1',
         R2 = 'LFQ intensity R2',
         R3 = 'LFQ intensity R3',
         R4 = 'LFQ intensity R4',
         R5 = 'LFQ intensity R5',
         C1 = 'LFQ intensity C1',
         C2 = 'LFQ intensity C2',
         C3 = 'LFQ intensity C3',
         C4 = 'LFQ intensity C4',
         C5 = 'LFQ intensity C5'
  )

protein_groups_iBAQ <- protein_groups %>%
  select(protein_ids = 'Majority protein IDs',
         R1 = 'iBAQ R1',
         R2 = 'iBAQ R2',
         R3 = 'iBAQ R3',
         R4 = 'iBAQ R4',
         R5 = 'iBAQ R5',
         C1 = 'iBAQ C1',
         C2 = 'iBAQ C2',
         C3 = 'iBAQ C3',
         C4 = 'iBAQ C4',
         C5 = 'iBAQ C5'
  )


## Check for negative binomial distribution

protein_groups_lfq %>%
  ggplot() +
  geom_histogram(aes(R1, fill = "R", color = "1"), bins = 100) +
  geom_histogram(aes(R2, fill = "R", color = "2"), bins = 100) +
  geom_histogram(aes(R3, fill = "R", color = "3"), bins = 100) +
  geom_histogram(aes(R4, fill = "R", color = "4"), bins = 100) +
  geom_histogram(aes(R5, fill = "R", color = "5"), bins = 100) +
  geom_histogram(aes(C1, fill = "C", color = "1"), bins = 100) +
  geom_histogram(aes(C2, fill = "C", color = "2"), bins = 100) +
  geom_histogram(aes(C3, fill = "C", color = "3"), bins = 100) +
  geom_histogram(aes(C4, fill = "C", color = "4"), bins = 100) +
  geom_histogram(aes(C5, fill = "C", color = "5"), bins = 100) +
  scale_y_continuous(trans = "log1p", breaks = c(0, 1, 1000, 2000, 4000, 6000))


#Create the initial **edgeR** object with all the values.

library(edgeR)
groups = c("R", "R", "R", "R", "R",
           "C", "C", "C", "C", "C") %>%
  factor(levels = c("R", "C")
  )
design = model.matrix(~0+groups)
contrasts = makeContrasts(
  H = groupsR - groupsC,
  levels = design
)
LFQ_edgeR <- DGEList(
  counts = protein_groups_lfq %>%
    mutate(protein_ids = make.unique(protein_ids)) %>%
    na.omit() %>%
    column_to_rownames("protein_ids"),
  group = c(1,1,1,1,1, 2,2,2,2,2)
)

#Filter out rows so that each row has at least 10 non-zero values.

LFQ_edgeR_10 <- LFQ_edgeR[ rowSums(cpm(LFQ_edgeR) == 0 ) <= (10 - 10) , ]
LFQ_edgeR_10 <- LFQ_edgeR_10 %>% calcNormFactors()
LFQ_edgeR_10 <- LFQ_edgeR_10 %>% estimateGLMCommonDisp(design = design)
## Following "decideTestsDGE(p.value = 0.10)" counts up, down, and unchanged genes (proteins) at 10% FDR
lrt_LFQ_10 <- glmLRT(glmFit(LFQ_edgeR_10, design), contrast = contrasts) %>% decideTestsDGE(p.value = 0.10) 
plotSmear(LFQ_edgeR_10, de.tags = rownames(lrt_LFQ_10)[lrt_LFQ_10 %>% as.logical()])

#Filter out rows so that each row has at least 8 non-zero values.

LFQ_edgeR_8 <- LFQ_edgeR[ rowSums(cpm(LFQ_edgeR) == 0 ) <= (10 - 8) , ]
LFQ_edgeR_8 <- LFQ_edgeR_8 %>% calcNormFactors()
LFQ_edgeR_8 <- LFQ_edgeR_8 %>% estimateGLMCommonDisp(design = design)
lrt_LFQ_8 <- glmLRT(glmFit(LFQ_edgeR_8, design), contrast = contrasts) %>% decideTestsDGE(p.value = 0.10)
plotSmear(LFQ_edgeR_8, de.tags = rownames(lrt_LFQ_8)[lrt_LFQ_8 %>% as.logical()])

#Filter out rows so that each row has at least 6 non-zero values.

LFQ_edgeR_6 <- LFQ_edgeR[ rowSums(cpm(LFQ_edgeR) == 0 ) <= (10 - 6) , ]
LFQ_edgeR_6 <- LFQ_edgeR_6 %>% calcNormFactors()
LFQ_edgeR_6 <- LFQ_edgeR_6 %>% estimateGLMCommonDisp(design = design)
lrt_LFQ_6 <- glmLRT(glmFit(LFQ_edgeR_6, design), contrast = contrasts) %>% decideTestsDGE(p.value = 0.10)
plotSmear(LFQ_edgeR_6, de.tags = rownames(lrt_LFQ_6)[lrt_LFQ_6 %>% as.logical()])

#Filter out rows so that each row has at least 4 non-zero values.

LFQ_edgeR_4 <- LFQ_edgeR[ rowSums(cpm(LFQ_edgeR) == 0 ) <= (10 - 4) , ]
LFQ_edgeR_4 <- LFQ_edgeR_4 %>% calcNormFactors()
LFQ_edgeR_4 <- LFQ_edgeR_4 %>% estimateGLMCommonDisp(design = design)
lrt_LFQ_4 <- glmLRT(glmFit(LFQ_edgeR_4, design), contrast = contrasts) %>% decideTestsDGE(p.value = 0.10)
plotSmear(LFQ_edgeR_4, de.tags = rownames(lrt_LFQ_4)[lrt_LFQ_4 %>% as.logical()])

#Filter out rows so that each row has at least 2 non-zero values.

LFQ_edgeR_2 <- LFQ_edgeR[ rowSums(cpm(LFQ_edgeR) == 0 ) <= (10 - 2) , ]
LFQ_edgeR_2 <- LFQ_edgeR_2 %>% calcNormFactors()
LFQ_edgeR_2 <- LFQ_edgeR_2 %>% estimateGLMCommonDisp(design = design)
lrt_LFQ_2 <- glmLRT(glmFit(LFQ_edgeR_2, design), contrast = contrasts) %>% decideTestsDGE(p.value = 0.10)
plotSmear(LFQ_edgeR_2, de.tags = rownames(lrt_LFQ_2)[lrt_LFQ_2 %>% as.logical()])

#Filter out rows so that each row has at least 1 non-zero values.

LFQ_edgeR_1 <- LFQ_edgeR[ rowSums(cpm(LFQ_edgeR) == 0 ) <= (10 - 1) , ]
LFQ_edgeR_1 <- LFQ_edgeR_1 %>% calcNormFactors()
LFQ_edgeR_1 <- LFQ_edgeR_1 %>% estimateGLMCommonDisp(design = design)
lrt_LFQ_1 <- glmLRT(glmFit(LFQ_edgeR_1, design), contrast = contrasts) %>% decideTestsDGE(p.value = 0.10)
plotSmear(LFQ_edgeR_1, de.tags = rownames(lrt_LFQ_1)[lrt_LFQ_1 %>% as.logical()])

## Make plots

make_pca_plot <- function(tp, tp_edger, ellipse = T) {
  tp %>%
    na.omit %>%
    filter(protein_ids %in% tp_edger) %>%
    mutate(protein_ids = abbreviate(protein_ids, minlength = 9)) %>%
    column_to_rownames("protein_ids") %>%
    as.matrix() %>%
    t %>%
    prcomp(center = T, scale = T) %>%
    ggbiplot( labels = c("R1", "R2", "R3", "R4", "R5","C1", "C2", "C3", "C4", "C5"),
              groups = c("R", "R", "R", "R", "R","C", "C", "C", "C", "C"),
              ellipse = ellipse) +
    theme_classic()
}
make_heatmap <- function(tp, tp_edger, ellipse = T) {
  tp %>%
    na.omit %>%
    filter(protein_ids %in% tp_edger) %>%
    mutate(protein_ids = abbreviate(protein_ids, minlength = 9)) %>%
    mutate_if(is.numeric, log1p) %>%
    column_to_rownames("protein_ids") %>%
    as.matrix() %>%
    pheatmap::pheatmap()
}
protein_groups_lfq %>%
  make_pca_plot(tp_edger = rownames(lrt_LFQ_10)[as.logical(lrt_LFQ_10)]) +
  ggtitle("At least 10 non-zero columns")
protein_groups_lfq %>%
  make_pca_plot(tp_edger = rownames(lrt_LFQ_8)[as.logical(lrt_LFQ_8)]) +
  ggtitle("At least 8 non-zero columns")
protein_groups_lfq %>%
  make_pca_plot(tp_edger = rownames(lrt_LFQ_6)[as.logical(lrt_LFQ_6)]) +
  ggtitle("At least 6 non-zero columns")
protein_groups_lfq %>%
  make_pca_plot(tp_edger = rownames(lrt_LFQ_4)[as.logical(lrt_LFQ_4)]) +
  ggtitle("At least 4 non-zero columns")
protein_groups_lfq %>%
  make_pca_plot(tp_edger = rownames(lrt_LFQ_2)[as.logical(lrt_LFQ_2)]) +
  ggtitle("At least 2 non-zero columns")
protein_groups_lfq %>%
  make_pca_plot(tp_edger = rownames(lrt_LFQ_1)[as.logical(lrt_LFQ_1)]) +
  ggtitle("At least 1 non-zero columns")



#Make a table for the column 10 or more.

table_LFQ_edgeR_10 <- (glmLRT(glmFit(LFQ_edgeR_10, design), contrast = contrasts) %>%
                        topTags(n=5000))$table %>%
  rownames_to_column("Majority protein IDs") %>%
  as_tibble() %>%
  merge(protein_groups %>%
          mutate(`Majority protein IDs` =
                   make.unique(`Majority protein IDs`))) %>%
  as_tibble() %>%
  arrange(FDR)

write.table(table_LFQ_edgeR_10 , file = "10x_full.csv", sep = ",", row.names = FALSE)


#Make a table for the column 8 or more.

table_LFQ_edgeR_8 <- (glmLRT(glmFit(LFQ_edgeR_8, design), contrast = contrasts) %>%
                         topTags(n=5000))$table %>%
  rownames_to_column("Majority protein IDs") %>%
  as_tibble() %>%
  merge(protein_groups %>%
          mutate(`Majority protein IDs` =
                   make.unique(`Majority protein IDs`))) %>%
  as_tibble() %>%
  arrange(FDR)

write.table(table_LFQ_edgeR_8 , file = "8x_full.csv", sep = ",", row.names = FALSE)


#Make a table for the column 6 or more.

table_LFQ_edgeR_6 <- (glmLRT(glmFit(LFQ_edgeR_6, design), contrast = contrasts) %>%
                        topTags(n=5000))$table %>%
  rownames_to_column("Majority protein IDs") %>%
  as_tibble() %>%
  merge(protein_groups %>%
          mutate(`Majority protein IDs` =
                   make.unique(`Majority protein IDs`))) %>%
  as_tibble() %>%
  arrange(FDR)

write.table(table_LFQ_edgeR_6 , file = "6x_full.csv", sep = ",", row.names = FALSE)


#Make a table for the column 4 or more.

table_LFQ_edgeR_4 <- (glmLRT(glmFit(LFQ_edgeR_4, design), contrast = contrasts) %>%
                        topTags(n=5000))$table %>%
  rownames_to_column("Majority protein IDs") %>%
  as_tibble() %>%
  merge(protein_groups %>%
          mutate(`Majority protein IDs` =
                   make.unique(`Majority protein IDs`))) %>%
  as_tibble() %>%
  arrange(FDR)

write.table(table_LFQ_edgeR_4 , file = "4x_full.csv", sep = ",", row.names = FALSE)



#Make a table for the column 2 or more.

table_LFQ_edgeR_2 <- (glmLRT(glmFit(LFQ_edgeR_2, design), contrast = contrasts) %>%
                        topTags(n=5000))$table %>%
  rownames_to_column("Majority protein IDs") %>%
  as_tibble() %>%
  merge(protein_groups %>%
          mutate(`Majority protein IDs` =
                   make.unique(`Majority protein IDs`))) %>%
  as_tibble() %>%
  arrange(FDR)

write.table(table_LFQ_edgeR_2 , file = "2x_full.csv", sep = ",", row.names = FALSE)


#Make a table for the column 1 or more.

table_LFQ_edgeR_1 <- (glmLRT(glmFit(LFQ_edgeR_1, design), contrast = contrasts) %>%
                        topTags(n=5000))$table %>%
  rownames_to_column("Majority protein IDs") %>%
  as_tibble() %>%
  merge(protein_groups %>%
          mutate(`Majority protein IDs` =
                   make.unique(`Majority protein IDs`))) %>%
  as_tibble() %>%
  arrange(FDR)

write.table(table_LFQ_edgeR_1 , file = "1x_full.csv", sep = ",", row.names = FALSE)




















