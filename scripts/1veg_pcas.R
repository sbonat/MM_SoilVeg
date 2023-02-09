library(ade4)
library(ape)
library(ecodist)
library(factoextra)
library(FactoMineR)
library(gclus)
library(ggforce)
library(ggplot2)
library(ggpubr)
library(ggsci)
library(ggrepel)
library(hrbrthemes)
library(lubridate)
library(readxl)
library(stats)
library(tidyverse)
library(vegan)
library(viridis)

# library(BiodiversityR) #makes R crash

#load data ####
veg0 <- read_excel("data/veg_clean.xlsx")
veg <- veg0 %>% 
       mutate(pre_or_post = gsub("Pre", "Baseline", pre_or_post),
              site = gsub("CR", "Site 1", site),
              site = gsub("MS", "Site 3", site),
              site = gsub("WP", "Site 2", site),
              treatment = gsub("HE", "Open", treatment),
              treatment = gsub("Nil", "Open", treatment),
              treatment = gsub("VE", "Vertebrate exclusion", treatment),
              treatment = gsub("IE", "Invertebrate exclusion", treatment),
              type = gsub("MM", "Mass mortality", type),
              type = gsub("SC", "Single carcass", type),
              type = gsub("C", "Control", type),
              position = transect, .after = pre_or_post,
              position = gsub("1", "Inside", position),
              position = gsub("2", "Inside", position),
              position = gsub("3", "Inside", position),
              position = gsub("4", "Outside", position),
              date = ymd(date),
             #Divide all cover variables by 100 so I have proportions
             deadmaterial = deadmaterial/100,
             bareground = bareground/100,
             cover = cover/100,
             forb = forb/100,
             shrub = shrub/100,
             mossfern = mossfern/100,
             grasses = grasses/100) %>% 
              filter(pre_or_post != "Year1") %>% 
            rename(debris = deadmaterial)
#View(veg)
# # Histograms and violin plots ####
# h1 <- hist(veg$debris) 
# h2 <- hist(veg$bareground)
# h3 <- hist(veg$vegcover)
# h4 <- hist(veg$forb)
# h5 <- hist(veg$grasssedgerush)
# h6 <- hist(veg$mossfern)
# h7 <- hist(veg$shrub)

#create violin plot for all above variables####

# #first, put data from wide to long format
# veg_long <- pivot_longer(veg, cols = c("debris", "bareground", "vegcover", "forb", "shrub", "mossfern", "grasssedgerush"), 
#                          names_to = "categories", values_to = "percentage")
# #Then select cover data
# veg_groundcover <- veg_long %>% 
#   filter(
#     categories != "forb", 
#     categories != "shrub", 
#     categories != "mossfern", 
#     categories != "grasssedgerush")
# #Then do violin plot for cover data
# #violin_coverdata <- ggplot(veg_groundcover, aes(x=categories, y= percentage, fill = treatment)) + 
# #  geom_violin()
# #violin_coverdata
# 
# #boxplot may be better
# groundcoverboxplot <- ggplot(veg_groundcover, aes(x=categories, y= percentage, fill = categories)) + 
#   geom_boxplot()
# groundcoverboxplot
# 
# 
# #Then select plant group data 
# veg_plantcat <- veg_long %>% 
#   filter(
#     categories != "bareground", 
#     categories != "vegcover",
#     categories != "debris")
# #And do a violin plot for the different plant categories
# plantboxplot <- ggplot(veg_plantcat, aes(x=categories, y= percentage, fill = categories)) + 
#   geom_boxplot()
# plantboxplot


# PCA1 all data####
pca1.data <- veg %>% 
             select(moisture, height, debris, bareground, cover, forb, shrub, mossfern,
                   grasses) 
#environmental variables
var1 <- veg %>% 
  select(id, site, type, treatment, pre_or_post, position)%>% 
    mutate(site = as.factor(site),
           type = as.factor(type),
           treatment = as.factor(treatment),
           position = as.factor(position),
           pre_or_post = as.factor(pre_or_post))
pca1 <- rda(pca1.data, scale = TRUE)
# summary(pca1)

# PCA1 eigenvalues and broken stick model to select PCA axes to consider ####
#     ev <- pca1$CA$eig #eigenvalues
#     ev[ev>mean(ev)] #apply Kaiser-Guttman criterion to select axes: only the axes above the mean computed should be interpreted
#     
#     #Broken stick model
#     n <- length(ev)
#     bsm <- data.frame(j=seq(1:n), p=0)
#     for (i in 2:n) {
#       bsm$p[i] = bsm$p[i-1] + (1/(n+1-i))
#       }
#     #Comparison of PCA result with broken stick model
#     par(mfrow = c(2,1))
#     barplot(ev, main = "Eigenvalues", col = "bisque", las=2)
#     abline(h=mean(ev), col = "red")
#     legend("topright", "Average eigenvalue", lwd = 1, col = 2, bty = "n")
#     barplot(t(cbind(100*ev/sum(ev), bsm$p[n:1])), beside = TRUE, main = "% variance", col=c("bisque",2), las=2)
# dev.off()

# PCA2 Subset data excluding baseline ####
pca2.data <- veg %>%
  filter(pre_or_post != "Baseline") %>%
  select(height, moisture, debris, bareground, cover, forb, shrub, mossfern,
         grasses)
var2 <- veg %>%
  filter(pre_or_post != "Baseline") %>%
  select(site, type, treatment, date, position)%>%
  mutate(site = as.factor(site),
         type = as.factor(type),
         treatment = as.factor(treatment),
         position = as.factor(position))
pca2 <- rda(pca2.data, scale = TRUE)
# summary(pca2)

# PCA2 eigenvalues and broken stick model to select PCA axes to consider ####
# ev <- pca2$CA$eig #eigenvalues
# ev[ev>mean(ev)] #apply Kaiser-Guttman criterion to select axes: only the axes above the mean computed should be interpreted
# 
# #Broken stick model
#     n <- length(ev)
#     bsm <- data.frame(j=seq(1:n), p=0)
#     for (i in 2:n) {
#       bsm$p[i] = bsm$p[i-1] + (1/(n+1-i))
#     }
#     #Comparison of PCA result with broken stick model
#     par(mfrow = c(2,1))
#     barplot(ev, main = "Eigenvalues", col = "bisque", las=2)
#     abline(h=mean(ev), col = "red")
#     legend("topright", "Average eigenvalue", lwd = 1, col = 2, bty = "n")
#     barplot(t(cbind(100*ev/sum(ev), bsm$p[n:1])), beside = TRUE, main = "% variance", col=c("bisque",2), las=2)
# #    dev.off()
# 

# PCA3 Subset data and eliminate outside plot transects####
pca3.data <- veg %>%
            filter(pre_or_post != "Baseline",
                   position != "Outside") %>%
            select(height, moisture, debris, bareground, cover, forb, shrub, mossfern,
                   grasses)

var3 <- veg %>%
            filter(pre_or_post != "Baseline",
                   position != "Outside") %>%
            select(site, type, treatment, date)%>%
            mutate(site = as.factor(site),
                   type = as.factor(type),
                   treatment = as.factor(treatment))
pca3 <- rda(pca3.data, scale = TRUE)
# summary(pca3, scaling = 1)

# PCA results ####
# results1 <- as.data.frame(pca1$CA$v[,c(1,2,3)])
# results1 <- results1 %>% mutate(veg = rownames(results1), .before = PC1,
#          pca_n = vector(mode = "character", length = 9),
#          pca_n = "pca1")
# results2 <- as.data.frame(pca2$CA$v[,c(1,2,3)])  
# results2 <- results2 %>%  mutate(veg = rownames(results2), .before = PC1,
#          pca_n = vector(mode = "character", length = 9),
#          pca_n = "pca2")
# results3 <- as.data.frame(pca3$CA$v[,c(1,2,3)])
# results3 <- results3 %>% mutate(veg = rownames(results3), .before = PC1,
#          pca_n = vector(mode = "character", length = 9),
#          pca_n = "pca3")
# 
# results_vegPCAs <- rbind(results1, results2, results3) #rbind is to stack the dfs on top of each other
# #column names must be the same
# Plots PCA1 ####
percent_explained1 <- format(round((100 * pca1$CA$eig/sum(pca1$CA$eig)), digits = 1), nsmall = 1, trim = TRUE)

colvec <- c("darkblue", "orange", "grey68")
attach(var1)

# Plots PCA1 Axis 1, 2  grouped by site####
png(filename = "figures/vegsoil/veg/pca1_bysite.png",
    width = 1415,
    height = 501,
    units = "px")

par(mfrow = c(1, 3))
ordiplot(pca1, choices = c(1,2), type = "n", 
         main = "Vegetation data grouped by site, Axis 1 and 2", 
         xlab = glue("PCA Axis 1 ({percent_explained1[1]}%)"),
         ylab = glue("PCA Axis 2 ({percent_explained1[2]}%)"),
         cex.lab=1.5, cex.axis=1.5, cex.main=1.5, cex.sub=1.5)
points(pca1, choices = c(1,2), col = colvec[site])
text(pca1, choices = c(1,2),display = "species", 
     scaling = 2, col = "red", cex = 1.5)
ordihull(pca1, choices = c(1,2), groups = site, col = colvec)
legend("topleft", legend = levels(site), bty = "n",
       col = colvec, pch = 21, pt.bg = colvec, cex = 1.5)
#PCA1 Axis 1, 3 grouped by site
ordiplot(pca1, choices = c(1,3), type = "n", 
         main = "Vegetation data grouped by site, Axis 1 and 3", 
         xlab = glue("PCA Axis 1 ({percent_explained1[1]}%)"),
         ylab = glue("PCA Axis 3 ({percent_explained1[3]}%)"),
         cex.lab=1.5, cex.axis=1.5, cex.main=1.5, cex.sub=1.5)
points(pca1, choices = c(1,3), col = colvec[site])
text(pca1, choices = c(1,3), display = "species", 
     scaling = 2, col = "red", cex = 1.5)
ordihull(pca1, choices = c(1,3), groups = site, col = colvec)
legend("topleft", legend = levels(site), bty = "n",
       col = colvec, pch = 21, pt.bg = colvec, cex = 1.5)

#PCA1 Axis 2, 3 grouped by site
ordiplot(pca1, choices = c(2,3), type = "n", 
         main = "Vegetation data grouped by site, Axis 2 and 3", 
         xlab = glue("PCA Axis 2 ({percent_explained1[2]}%)"),
         ylab = glue("PCA Axis 3 ({percent_explained1[3]}%)"),
         cex.lab=1.5, cex.axis=1.5, cex.main=1.5, cex.sub=1.5)
points(pca1, choices = c(2,3), col = colvec[site])
text(pca1, choices = c(2,3),display = "species", 
     scaling = 2, col = "red", cex = 1.5)
ordihull(pca1, choices = c(2,3), groups = site, col = colvec)
legend("topleft", legend = levels(site), bty = "n",
       col = colvec, pch = 21, pt.bg = colvec, cex = 1.5)
dev.off()

# Plots PCA1 Axis 1, 2 grouped by type####    
png(filename = "figures/vegsoil/veg/pca1_bytype.png",
    width = 1415,
    height = 501,
    units = "px")

par(mfrow = c(1, 3))

ordiplot(pca1, choices = c(1,2), type = "n", 
         main = "Vegetation data grouped by type, Axis 1 and 2", 
         xlab = glue("PCA Axis 1 ({percent_explained1[1]}%)"),
         ylab = glue("PCA Axis 2 ({percent_explained1[2]}%)"),
         cex.lab=1.5, cex.axis=1.5, cex.main=1.5, cex.sub=1.5)
points(pca1, choices = c(1,2), col = colvec[type])
text(pca1, choices = c(1,2),display = "species", 
     scaling = 2, col = "red", cex = 1.5)
ordihull(pca1, choices = c(1,2), groups = type, col = colvec)
legend("topleft", legend = levels(type), bty = "n",
       col = colvec, pch = 21, pt.bg = colvec, cex = 1.5)
#PCA1 type Axis 1, 3 grouped by type
ordiplot(pca1, choices = c(1,3), type = "n", 
         main = "Vegetation data grouped by type, Axis 1 and 3", 
         xlab = glue("PCA Axis 1 ({percent_explained1[1]}%)"),
         ylab = glue("PCA Axis 3 ({percent_explained1[3]}%)"),
         cex.lab=1.5, cex.axis=1.5, cex.main=1.5, cex.sub=1.5)
points(pca1, choices = c(1,3), col = colvec[type])
text(pca1, choices = c(1,3),display = "species", 
     scaling = 2, col = "red", cex = 1.5)
ordihull(pca1, choices = c(1,3), groups = type, col = colvec)
legend("topleft", legend = levels(type), bty = "n",
       col = colvec, pch = 21, pt.bg = colvec, cex = 1.5)

#PCA1 type Axis 2, 3 grouped by type
ordiplot(pca1, choices = c(2,3), type = "n", 
         main = "Vegetation data grouped by type, Axis 2 and 3", 
         xlab = glue("PCA Axis 2 ({percent_explained1[2]}%)"),
         ylab = glue("PCA Axis 3 ({percent_explained1[3]}%)"),
         cex.lab=1.5, cex.axis=1.5, cex.main=1.5, cex.sub=1.5)
points(pca1, choices = c(2,3), col = colvec[type])
text(pca1, choices = c(2,3),display = "species", 
     scaling = 2, col = "red", cex = 1.5)
ordihull(pca1, choices = c(2,3), groups = type, col = colvec)
legend("topleft", legend = levels(type), bty = "n",
       col = colvec, pch = 21, pt.bg = colvec, cex = 1.5)

dev.off()
# Plots PCA1 Axis 1, 2 grouped by treatment####
png(filename = "figures/vegsoil/veg/pca1_bytreat.png",
    width = 1415,
    height = 501,
    units = "px")

par(mfrow = c(1, 3))
ordiplot(pca1, choices = c(1,2), type = "n", 
         main = "Vegetation data grouped by treatment, Axis 1 and 2", 
         xlab = glue("PCA Axis 1 ({percent_explained1[1]}%)"),
         ylab = glue("PCA Axis 2 ({percent_explained1[2]}%)"),
         cex.lab=1.5, cex.axis=1.5, cex.main=1.5, cex.sub=1.5)
points(pca1, choices = c(1,2), col = colvec[treatment])
text(pca1, choices = c(1,2),display = "species", 
     scaling = 2, col = "red", cex = 1.5)
ordihull(pca1, choices = c(1,2), groups = treatment, col = colvec)
legend("topleft", legend = levels(treatment), bty = "n",
       col = colvec, pch = 21, pt.bg = colvec, cex = 1.5)
#PCA1 type Axis 1, 3 grouped by treatment
ordiplot(pca1, choices = c(1,3), type = "n", 
         main = "Vegetation data grouped by treatment, Axis 1 and 3", 
         xlab = glue("PCA Axis 1 ({percent_explained1[1]}%)"),
         ylab = glue("PCA Axis 3 ({percent_explained1[3]}%)"),
         cex.lab=1.5, cex.axis=1.5, cex.main=1.5, cex.sub=1.5)
points(pca1, choices = c(1,3), col = colvec[treatment])
text(pca1, choices = c(1,3),display = "species", 
     scaling = 2, col = "red", cex = 1.5)
ordihull(pca1, choices = c(1,3), groups = treatment, col = colvec)
legend("topleft", legend = levels(treatment), bty = "n",
       col = colvec, pch = 21, pt.bg = colvec, cex = 1.5)

#PCA1 type Axis 2, 3 grouped by treatment
ordiplot(pca1, choices = c(2,3), type = "n", 
         main = "Vegetation data grouped by treatment, Axis 2 and 3", 
         xlab = glue("PCA Axis 2 ({percent_explained1[2]}%)"),
         ylab = glue("PCA Axis 3 ({percent_explained1[3]}%)"),
         cex.lab=1.5, cex.axis=1.5, cex.main=1.5, cex.sub=1.5)
points(pca1, choices = c(2,3), col = colvec[treatment])
text(pca1, choices = c(2,3),display = "species", 
     scaling = 2, col = "red", cex = 1.5)
ordihull(pca1, choices = c(2,3), groups = treatment, col = colvec)
legend("topleft", legend = levels(treatment), bty = "n",
       col = colvec, pch = 21, pt.bg = colvec, cex = 1.5)
dev.off()
# Plots PCA1 Axis 1, 2 grouped by position####
png(filename = "figures/vegsoil/veg/pca1_byposition.png",
    width = 1415,
    height = 501,
    units = "px")

par(mfrow = c(1, 3))
ordiplot(pca1, choices = c(1,2), type = "n", 
         main = "Vegetation data inside vs outside plots, Axis 1 and 2", 
         xlab = glue("PCA Axis 1 ({percent_explained1[1]}%)"),
         ylab = glue("PCA Axis 2 ({percent_explained1[2]}%)"),
         cex.lab=1.5, cex.axis=1.5, cex.main=1.5, cex.sub=1.5)
points(pca1, choices = c(1,2), col = colvec[position])
text(pca1, choices = c(1,2),display = "species", 
     scaling = 2, col = "red", cex = 1.5)
ordihull(pca1, choices = c(1,2), groups = position, col = colvec)
legend("topleft", legend = levels(position), bty = "n",
       col = colvec, pch = 21, pt.bg = colvec, cex = 1.5)
#PCA1 type Axis 1, 3 grouped by position
ordiplot(pca1, choices = c(1,3), type = "n", 
         main = "Vegetation data inside vs outside plots, Axis 1 and 3", 
         xlab = glue("PCA Axis 1 ({percent_explained1[1]}%)"),
         ylab = glue("PCA Axis 3 ({percent_explained1[3]}%)"),
         cex.lab=1.5, cex.axis=1.5, cex.main=1.5, cex.sub=1.5)
points(pca1, choices = c(1,3), col = colvec[position])
text(pca1, choices = c(1,3),display = "species", 
     scaling = 2, col = "red", cex = 1.5)
ordihull(pca1, choices = c(1,3), groups = position, col = colvec)
legend("topleft", legend = levels(position), bty = "n",
       col = colvec, pch = 21, pt.bg = colvec, cex = 1.5)

#PCA1 type Axis 2, 3 grouped by position
ordiplot(pca1, choices = c(2,3), type = "n", 
         main = "Vegetation data inside vs outside plots, Axis 2 and 3", 
         xlab = glue("PCA Axis 2 ({percent_explained1[2]}%)"),
         ylab = glue("PCA Axis 3 ({percent_explained1[3]}%)"),
         cex.lab=1.5, cex.axis=1.5, cex.main=1.5, cex.sub=1.5)
points(pca1, choices = c(2,3), col = colvec[position])
text(pca1, choices = c(2,3),display = "species", 
     scaling = 2, col = "red", cex = 1.5)
ordihull(pca1, choices = c(2,3), groups = position, col = colvec)
legend("topleft", legend = levels(position), bty = "n",
       col = colvec, pch = 21, pt.bg = colvec, cex = 1.5)
dev.off()
# Plots PCA1 Axis 1, 2 grouped by pre_or_post####
png(filename = "figures/vegsoil/veg/pca1_preorpost.png",
    width = 1415,
    height = 501,
    units = "px")
par(mfrow = c(1, 3))
ordiplot(pca1, choices = c(1,2), type = "n", 
         main = "Vegetation data: baseline vs post data, Axis 1 and 2", 
         xlab = glue("PCA Axis 1 ({percent_explained1[1]}%)"),
         ylab = glue("PCA Axis 2 ({percent_explained1[2]}%)"),
         cex.lab=1.5, cex.axis=1.5, cex.main=1.5, cex.sub=1.5)
points(pca1, choices = c(1,2), col = colvec[pre_or_post])
text(pca1, choices = c(1,2),display = "species", 
     scaling = 2, col = "red", cex = 1.5)
ordihull(pca1, choices = c(1,2), groups = pre_or_post, col = colvec)
legend("topleft", legend = levels(pre_or_post), bty = "n",
       col = colvec, pch = 21, pt.bg = colvec, cex = 1.5)
#PCA1 type Axis 1, 3 grouped by pre_or_post
ordiplot(pca1, choices = c(1,3), type = "n", 
         main = "Vegetation data: baseline vs post data, Axis 1 and 3", 
         xlab = glue("PCA Axis 1 ({percent_explained1[1]}%)"),
         ylab = glue("PCA Axis 3 ({percent_explained1[3]}%)"),
         cex.lab=1.5, cex.axis=1.5, cex.main=1.5, cex.sub=1.5)
points(pca1, choices = c(1,3), col = colvec[pre_or_post])
text(pca1, choices = c(1,3),display = "species", 
     scaling = 2, col = "red", cex = 1.5)
ordihull(pca1, choices = c(1,3), groups = pre_or_post, col = colvec)
legend("topleft", legend = levels(pre_or_post), bty = "n",
       col = colvec, pch = 21, pt.bg = colvec, cex = 1.5)

#PCA1 type Axis 2, 3 grouped by pre_or_post
ordiplot(pca1, choices = c(2,3), type = "n", 
         main = "Vegetation data: baseline vs post data, Axis 2 and 3", 
         xlab = glue("PCA Axis 2 ({percent_explained1[2]}%)"),
         ylab = glue("PCA Axis 3 ({percent_explained1[3]}%)"),
         cex.lab=1.5, cex.axis=1.5, cex.main=1.5, cex.sub=1.5)
points(pca1, choices = c(2,3), col = colvec[pre_or_post])
text(pca1, choices = c(2,3),display = "species", 
     scaling = 2, col = "red", cex = 1.5)
ordihull(pca1, choices = c(2,3), groups = pre_or_post, col = colvec)
legend("topleft", legend = levels(pre_or_post), bty = "n",
       col = colvec, pch = 21, pt.bg = colvec, cex = 1.5)
dev.off()  
# Plots for PCA2 ####
percent_explained2 <- format(round((100 * pca2$CA$eig/sum(pca2$CA$eig)), digits = 1), nsmall = 1, trim = TRUE)
colvec <- c("darkblue", "orange", "grey68")

detach(var1)
attach(var2)
# Plots pca2 Axis 1, 2 grouped by site####
png(filename = "figures/vegsoil/veg/pca2_bysite.png",
    width = 1415,
    height = 501,
    units = "px")

par(mfrow = c(1, 3))
ordiplot(pca2, choices = c(1,2), type = "n", 
         main = "Vegetation data grouped by site, Axis 1 and 2", 
         xlab = glue("PCA Axis 1 ({percent_explained2[1]}%)"),
         ylab = glue("PCA Axis 2 ({percent_explained2[2]}%)"),
         cex.lab=1.5, cex.axis=1.5, cex.main=1.5, cex.sub=1.5)
points(pca2, choices = c(1,2), col = colvec[site])
text(pca2, choices = c(1,2),display = "species", 
     scaling = 2, col = "red", cex = 1.5)
ordihull(pca2, choices = c(1,2), groups = site, col = colvec)
legend("topleft", legend = levels(site), bty = "n",
       col = colvec, pch = 21, pt.bg = colvec, cex = 1.5)
# Plots pca2 Axis 1, 3 grouped by site
ordiplot(pca2, choices = c(1,3), type = "n", 
         main = "Vegetation data grouped by site, Axis 1 and 3", 
         xlab = glue("PCA Axis 1 ({percent_explained2[1]}%)"),
         ylab = glue("PCA Axis 3 ({percent_explained2[3]}%)"),
         cex.lab=1.5, cex.axis=1.5, cex.main=1.5, cex.sub=1.5)
points(pca2, choices = c(1,3), col = colvec[site])
text(pca2, choices = c(1,3),display = "species", 
     scaling = 2, col = "red", cex = 1.5)
ordihull(pca2, choices = c(1,3), groups = site, col = colvec)
legend("topleft", legend = levels(site), bty = "n",
       col = colvec, pch = 21, pt.bg = colvec, cex = 1.5)

# Plots pca2 Axis 2, 3 grouped by site
ordiplot(pca2, choices = c(2,3), type = "n", 
         main = "Vegetation data grouped by site, Axis 2 and 3", 
         xlab = glue("PCA Axis 2 ({percent_explained2[2]}%)"),
         ylab = glue("PCA Axis 3 ({percent_explained2[3]}%)"),
         cex.lab=1.5, cex.axis=1.5, cex.main=1.5, cex.sub=1.5)
points(pca2, choices = c(2,3), col = colvec[site])
text(pca2, choices = c(2,3),display = "species", 
     scaling = 2, col = "red", cex = 1.5)
ordihull(pca2, choices = c(2,3), groups = site, col = colvec)
legend("topleft", legend = levels(site), bty = "n",
       col = colvec, pch = 21, pt.bg = colvec, cex = 1.5)
dev.off()

# Plots pca2 Axis 1, 2 grouped by type####   
png(filename = "figures/vegsoil/veg/pca2_bytype.png",
    width = 1415,
    height = 501,
    units = "px")
par(mfrow = c(1, 3))
ordiplot(pca2, choices = c(1,2), type = "n", 
         main = "Vegetation data grouped by type, Axis 1 and 2", 
         xlab = glue("PCA Axis 1 ({percent_explained2[1]}%)"),
         ylab = glue("PCA Axis 2 ({percent_explained2[2]}%)"),
         cex.lab=1.5, cex.axis=1.5, cex.main=1.5, cex.sub=1.5)
points(pca2, choices = c(1,2), col = colvec[type])
text(pca2, choices = c(1,2),display = "species", 
     scaling = 2, col = "red", cex = 1.5)
ordihull(pca2, choices = c(1,2), groups = type, col = colvec)
legend("topleft", legend = levels(type), bty = "n",
       col = colvec, pch = 21, pt.bg = colvec, cex = 1.5)
# Plots pca2 Axis 1, 3 grouped by type
ordiplot(pca2, choices = c(1,3), type = "n", 
         main = "Vegetation data grouped by type, Axis 1 and 3", 
         xlab = glue("PCA Axis 1 ({percent_explained2[1]}%)"),
         ylab = glue("PCA Axis 3 ({percent_explained2[3]}%)"),
         cex.lab=1.5, cex.axis=1.5, cex.main=1.5, cex.sub=1.5)
points(pca2, choices = c(1,3), col = colvec[type])
text(pca2, choices = c(1,3),display = "species", 
     scaling = 2, col = "red", cex = 1.5)
ordihull(pca2, choices = c(1,3), groups = type, col = colvec)
legend("topleft", legend = levels(type), bty = "n",
       col = colvec, pch = 21, pt.bg = colvec, cex = 1.5)
# Plots pca2 Axis 2, 3 grouped by type
ordiplot(pca2, choices = c(2,3), type = "n", 
         main = "Vegetation data grouped by type, Axis 2 and 3", 
         xlab = glue("PCA Axis 2 ({percent_explained2[2]}%)"),
         ylab = glue("PCA Axis 3 ({percent_explained2[3]}%)"),
         cex.lab=1.5, cex.axis=1.5, cex.main=1.5, cex.sub=1.5)
points(pca2, choices = c(2,3), col = colvec[type])
text(pca2, choices = c(2,3),display = "species", 
     scaling = 2, col = "red", cex = 1.5)
ordihull(pca2, choices = c(2,3), groups = type, col = colvec)
legend("bottomleft", legend = levels(type), bty = "n",
       col = colvec, pch = 21, pt.bg = colvec, cex = 1.5)
dev.off()    
# Plots pca2 Axis 1, 2 grouped by treatment####
png(filename = "figures/vegsoil/veg/pca2_bytreat.png",
    width = 1415,
    height = 501,
    units = "px")

par(mfrow = c(1, 3))
ordiplot(pca2, choices = c(1,2), type = "n", 
         main = "Vegetation data grouped by treatment, Axis 1 and 2", 
         xlab = glue("PCA Axis 1 ({percent_explained2[1]}%)"),
         ylab = glue("PCA Axis 2 ({percent_explained2[2]}%)"),
         cex.lab=1.5, cex.axis=1.5, cex.main=1.5, cex.sub=1.5)
points(pca2, choices = c(1,2), col = colvec[treatment])
text(pca2, choices = c(1,2),display = "species", 
     scaling = 2, col = "red", cex = 1.5)
ordihull(pca2, choices = c(1,2), groups = treatment, col = colvec)
legend("topleft", legend = levels(treatment), bty = "n",
       col = colvec, pch = 21, pt.bg = colvec, cex = 1.5)
# Plots pca2 Axi, 3 grouped by treatment
ordiplot(pca2, choices = c(1,3), type = "n", 
         main = "Vegetation data grouped by treatment, Axis 1 and 3", 
         xlab = glue("PCA Axis 1 ({percent_explained2[1]}%)"),
         ylab = glue("PCA Axis 3 ({percent_explained2[3]}%)"),
         cex.lab=1.5, cex.axis=1.5, cex.main=1.5, cex.sub=1.5)
points(pca2, choices = c(1,3), col = colvec[treatment])
text(pca2, choices = c(1,3),display = "species", 
     scaling = 2, col = "red", cex = 1.5)
ordihull(pca2, choices = c(1,3), groups = treatment, col = colvec)
legend("topleft", legend = levels(treatment), bty = "n",
       col = colvec, pch = 21, pt.bg = colvec, cex = 1.5)

# Plots pca2 Axi, 3 grouped by treatment
ordiplot(pca2, choices = c(2,3), type = "n", 
         main = "Vegetation data grouped by treatment, Axis 2 and 3", 
         xlab = glue("PCA Axis 2 ({percent_explained2[2]}%)"),
         ylab = glue("PCA Axis 3 ({percent_explained2[3]}%)"),
         cex.lab=1.5, cex.axis=1.5, cex.main=1.5, cex.sub=1.5)
points(pca2, choices = c(2,3), col = colvec[treatment])
text(pca2, choices = c(2,3),display = "species", 
     scaling = 2, col = "red", cex = 1.5)
ordihull(pca2, choices = c(2,3), groups = treatment, col = colvec)
legend("bottomleft", legend = levels(treatment), bty = "n",
       col = colvec, pch = 21, pt.bg = colvec, cex = 1.5)
dev.off()
# Plots pca2 Axis 1, 2 grouped by position####
png(filename = "figures/vegsoil/veg/pca2_byposition.png",
    width = 1415,
    height = 501,
    units = "px")

par(mfrow = c(1, 3))
ordiplot(pca2, choices = c(1,2), type = "n", 
         main = "Soil nutrients grouped by position, Axis 1 and 2", 
         xlab = glue("PCA Axis 1 ({percent_explained2[1]}%)"),
         ylab = glue("PCA Axis 2 ({percent_explained2[2]}%)"),
         cex.lab=1.5, cex.axis=1.5, cex.main=1.5, cex.sub=1.5)
points(pca2, choices = c(1,2), col = colvec[position])
text(pca2, choices = c(1,2),display = "species", 
     scaling = 2, col = "red", cex = 1.5)
ordihull(pca2, choices = c(1,2), groups = position, col = colvec)
legend("topleft", legend = levels(position), bty = "n",
       col = colvec, pch = 21, pt.bg = colvec, cex = 1.5)    
# Plots pca2 Axis 1, 3 grouped by position
ordiplot(pca2, choices = c(1,3), type = "n", 
         main = "Soil nutrients grouped by position, Axis 1 and 3", 
         xlab = glue("PCA Axis 1 ({percent_explained2[1]}%)"),
         ylab = glue("PCA Axis 3 ({percent_explained2[3]}%)"),
         cex.lab=1.5, cex.axis=1.5, cex.main=1.5, cex.sub=1.5)
points(pca2, choices = c(1,3), col = colvec[position])
text(pca2, choices = c(1,3),display = "species", 
     scaling = 2, col = "red", cex = 1.5)
ordihull(pca2, choices = c(1,3), groups = position, col = colvec)
legend("topleft", legend = levels(position), bty = "n",
       col = colvec, pch = 21, pt.bg = colvec, cex = 1.5)
# Plots pca2 Axis 2, 3 grouped by position
ordiplot(pca2, choices = c(2,3), type = "n", 
         main = "Vegetation data inside vs outside plots, Axis 2 and 3", 
         xlab = glue("PCA Axis 2 ({percent_explained2[2]}%)"),
         ylab = glue("PCA Axis 3 ({percent_explained2[3]}%)"),
         cex.lab=1.5, cex.axis=1.5, cex.main=1.5, cex.sub=1.5)
points(pca2, choices = c(2,3), col = colvec[position])
text(pca2, choices = c(2,3),display = "species", 
     scaling = 2, col = "red", cex = 1.5)
ordihull(pca2, choices = c(2,3), groups = position, col = colvec)
legend("topleft", legend = levels(position), bty = "n",
       col = colvec, pch = 21, pt.bg = colvec, cex = 1.5)
dev.off()
# PCA3 plots ####    
percent_explained3<- format(round((100 * pca3$CA$eig/sum(pca3$CA$eig)), digits = 1), nsmall = 1, trim = TRUE)
colvec <- c("darkblue", "orange", "grey68")

detach(var2)
attach(var3)
# Plots pca3 Axis 1, 2 grouped by site####
png(filename = "figures/vegsoil/veg/pca3_bysite.png",
    width = 1415,
    height = 501,
    units = "px")

par(mfrow = c(1, 3))
ordiplot(pca3, choices = c(1,2), type = "n", 
         main = "Vegetation data grouped by site, Axis 1 and 2", 
         xlab = glue("PCA Axis 1 ({percent_explained3[1]}%)"),
         ylab = glue("PCA Axis 2 ({percent_explained3[2]}%)"),
         cex.lab=1.5, cex.axis=1.5, cex.main=1.5, cex.sub=1.5)
points(pca3, choices = c(1,2), col = colvec[site])
text(pca3, choices = c(1,2),display = "species", 
     scaling = 2, col = "red", cex = 1.5)
ordihull(pca3, choices = c(1,2), groups = site, col = colvec)
legend("topright", legend = levels(site), bty = "n",
       col = colvec, pch = 21, pt.bg = colvec, cex = 1.5)
# Plots pca3 Axis 1, 3 grouped by site
ordiplot(pca3, choices = c(1,3), type = "n", 
         main = "Vegetation data grouped by site, Axis 1 and 3", 
         xlab = glue("PCA Axis 1 ({percent_explained3[1]}%)"),
         ylab = glue("PCA Axis 3 ({percent_explained3[3]}%)"),
         cex.lab=1.5, cex.axis=1.5, cex.main=1.5, cex.sub=1.5)
points(pca3, choices = c(1,3), col = colvec[site])
text(pca3, choices = c(1,3),display = "species", 
     scaling = 2, col = "red", cex = 1.5)
ordihull(pca3, choices = c(1,3), groups = site, col = colvec)
legend("topright", legend = levels(site), bty = "n",
       col = colvec, pch = 21, pt.bg = colvec, cex = 1.5)

# Plots pca3 Axis 2, 3 grouped by site
ordiplot(pca3, choices = c(2,3), type = "n", 
         main = "Vegetation data grouped by site, Axis 2 and 3", 
         xlab = glue("PCA Axis 2 ({percent_explained3[2]}%)"),
         ylab = glue("PCA Axis 3 ({percent_explained3[3]}%)"),
         cex.lab=1.5, cex.axis=1.5, cex.main=1.5, cex.sub=1.5)
points(pca3, choices = c(2,3), col = colvec[site])
text(pca3, choices = c(2,3),display = "species", 
     scaling = 2, col = "red", cex = 1.5)
ordihull(pca3, choices = c(2,3), groups = site, col = colvec)
legend("topright", legend = levels(site), bty = "n",
       col = colvec, pch = 21, pt.bg = colvec, cex = 1.5)
dev.off()  
# Plots pca3 Axis 1, 2 grouped by type####  
png(filename = "figures/vegsoil/veg/pca3_bytype.png",
    width = 1415,
    height = 501,
    units = "px")

par(mfrow = c(1, 3))

ordiplot(pca3, choices = c(1,2), type = "n", 
         main = "Vegetation data grouped by type, Axis 1 and 2", 
         xlab = glue("PCA Axis 1 ({percent_explained3[1]}%)"),
         ylab = glue("PCA Axis 2 ({percent_explained3[2]}%)"),
         cex.lab=1.5, cex.axis=1.5, cex.main=1.5, cex.sub=1.5)
points(pca3, choices = c(1,2), col = colvec[type])
text(pca3, choices = c(1,2),display = "species", 
     scaling = 2, col = "red", cex = 1.5)
ordihull(pca3, choices = c(1,2), groups = type, col = colvec)
legend("topright", legend = levels(type), bty = "n",
       col = colvec, pch = 21, pt.bg = colvec, cex = 1.5)
# Plots pca3 Axis 1, 3 grouped by type
ordiplot(pca3, choices = c(1,3), type = "n", 
         main = "Vegetation data grouped by type, Axis 1 and 3", 
         xlab = glue("PCA Axis 1 ({percent_explained3[1]}%)"),
         ylab = glue("PCA Axis 3 ({percent_explained3[3]}%)"),
         cex.lab=1.5, cex.axis=1.5, cex.main=1.5, cex.sub=1.5)
points(pca3, choices = c(1,3), col = colvec[type])
text(pca3, choices = c(1,3),display = "species", 
     scaling = 2, col = "red", cex = 1.5)
ordihull(pca3, choices = c(1,3), groups = type, col = colvec)
legend("topright", legend = levels(type), bty = "n",
       col = colvec, pch = 21, pt.bg = colvec, cex = 1.5)

# Plots pca3 Axis 2, 3 grouped by type
ordiplot(pca3, choices = c(2,3), type = "n", 
         main = "Vegetation data grouped by type, Axis 2 and 3", 
         xlab = glue("PCA Axis 2 ({percent_explained3[2]}%)"),
         ylab = glue("PCA Axis 3 ({percent_explained3[3]}%)"),
         cex.lab=1.5, cex.axis=1.5, cex.main=1.5, cex.sub=1.5)
points(pca3, choices = c(2,3), col = colvec[type])
text(pca3, choices = c(2,3),display = "species", 
     scaling = 2, col = "red", cex = 1.5)
ordihull(pca3, choices = c(2,3), groups = type, col = colvec)
legend("topright", legend = levels(type), bty = "n",
       col = colvec, pch = 21, pt.bg = colvec, cex = 1.5)
dev.off()   
# Plots pca3 Axis 1, 2 grouped by treatment####
png(filename = "figures/vegsoil/veg/pca3_bytreat.png",
    width = 1415,
    height = 501,
    units = "px")

par(mfrow = c(1, 3))

ordiplot(pca3, choices = c(1,2), type = "n", 
         main = "Vegetation data grouped by treatment, Axis 1 and 2", 
         xlab = glue("PCA Axis 1 ({percent_explained3[1]}%)"),
         ylab = glue("PCA Axis 2 ({percent_explained3[2]}%)"),
         cex.lab=1.5, cex.axis=1.5, cex.main=1.5, cex.sub=1.5)
points(pca3, choices = c(1,2), col = colvec[treatment])
text(pca3, choices = c(1,2),display = "species", 
     scaling = 2, col = "red", cex = 1.5)
ordihull(pca3, choices = c(1,2), groups = treatment, col = colvec)
legend("topright", legend = levels(treatment), bty = "n",
       col = colvec, pch = 21, pt.bg = colvec, cex = 1.5)
# Plots pca3 Axis 1, 3 grouped by treatment
ordiplot(pca3, choices = c(1,3), type = "n", 
         main = "Vegetation data grouped by treatment, Axis 1 and 3", 
         xlab = glue("PCA Axis 1 ({percent_explained3[1]}%)"),
         ylab = glue("PCA Axis 3 ({percent_explained3[3]}%)"),
         cex.lab=1.5, cex.axis=1.5, cex.main=1.5, cex.sub=1.5)
points(pca3, choices = c(1,3), col = colvec[treatment])
text(pca3, choices = c(1,3),display = "species", 
     scaling = 2, col = "red", cex = 1.5)
ordihull(pca3, choices = c(1,3), groups = treatment, col = colvec)
legend("topright", legend = levels(treatment), bty = "n",
       col = colvec, pch = 21, pt.bg = colvec, cex = 1.5)

# Plots pca3 Axis 2, 3 grouped by treatment
ordiplot(pca3, choices = c(2,3), type = "n", 
         main = "Vegetation data grouped by treatment, Axis 2 and 3", 
         xlab = glue("PCA Axis 2 ({percent_explained3[2]}%)"),
         ylab = glue("PCA Axis 3 ({percent_explained3[3]}%)"),
         cex.lab=1.5, cex.axis=1.5, cex.main=1.5, cex.sub=1.5)
points(pca3, choices = c(2,3), col = colvec[treatment])
text(pca3, choices = c(2,3),display = "species", 
     scaling = 2, col = "red", cex = 1.5)
ordihull(pca3, choices = c(2,3), groups = treatment, col = colvec)
legend("topright", legend = levels(treatment), bty = "n",
       col = colvec, pch = 21, pt.bg = colvec, cex = 1.5)
dev.off()