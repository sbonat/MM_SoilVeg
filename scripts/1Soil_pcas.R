#packages to load####
library(ade4)
library(ape)
library(cowplot)
library(ecodist)
library(ggplot2)
library(ggpubr)
library(glue)
library(lubridate)
library(readxl)
library(tidytidbits)
library(tidyverse)
library(vegan)
library(viridis)
library(concaveman)
library(pvclust)
library(RColorBrewer)

# Load data####
# setwd("C:/Users/sbon6737/Desktop/data/MME")
soil0 <- read_excel("data/soil_probes_clean.xlsx")
soil <- soil0 %>% 
        mutate(treatment = gsub("HE", "Open", treatment),
               treatment = gsub("Nil", "Open", treatment),
               treatment = gsub("VE", "Vertebrate exclusion", treatment),
               treatment = gsub("IE", "Invertebrate exclusion", treatment),
               type = gsub("MM", "Mass mortality", type),
               type = gsub("SC", "Single carcass", type),
               site = gsub("CR", "Site 1", site),
               site = gsub("MS", "Site 3", site),
               site = gsub("WP", "Site 2", site),
               pre_or_post = gsub("Pre", "Baseline", pre_or_post),
               B = NULL, Cd = NULL, Cu = NULL, Pb = NULL, Zn = NULL, Al = NULL, S = NULL,
               '#anion' = NULL, '#cation' = NULL, burial_date = NULL, retrieval_date = NULL) %>% 
        filter(pre_or_post != "Yr1") #selecting out Yr1 results as not relevant

# View(soil)

# Exploratory plots ####
# hist(soil$NO3_N)
# hist(soil$NH4_N)
# hist(soil$Ca)
# hist(soil$Mg)
# hist(soil$K)
# hist(soil$P)
# hist(soil$Fe)
# hist(soil$Mn)
# hist(soil$Cu)
# hist(soil$Zn)
# hist(soil$B)
# hist(soil$S)
# hist(soil$Pb)
# hist(soil$Al)
# hist(soil$Cd)
# 
# long_soil <- pivot_longer(soil, cols = c("NO3_N","NH4_N", "Ca","Mg", "K", "P" ,"Fe", "Mn" , "Cu" , "Zn", "B"  , "S", "Pb", "Al", "Cd"),
#                           names_to = "Ion", values_to = "Value")
# long_soil
# 
# boxplot <- ggplot(long_soil, aes(x=Ion, y= Value, fill = Ion))+
#             geom_boxplot()
# boxplot
# 
# group1 <- long_soil %>% 
#           filter(Ion != "Ca", Ion != "K", Ion != "Mg", Ion != "NH4_N",
#                  Ion != "Mn", Ion != "Fe", Ion != "NO3_N", 
#                  Ion !=   "Al", Ion !=  "P", Ion !=   "S",Ion != "Zn")
# group1
# 
# group2 <- long_soil %>% 
#   filter(Ion != "Mn", Ion != "Fe", Ion != "NO3_N",
#          Ion !=   "Al", Ion !=  "P", Ion !=   "S",Ion != "Zn",
#          Ion !="B", Ion !=  "Cd", Ion !=   "Cu",Ion != "Pb"
#          )
# group2
# 
# group3 <- long_soil %>% 
#   filter(Ion != "Ca", Ion != "K", Ion != "Mg", Ion != "NH4_N", 
#          Ion !=   "Al", Ion !=  "P", Ion !=   "S",Ion != "Zn",
#          Ion !="B", Ion !=  "Cd", Ion !=   "Cu",Ion != "Pb")
# group3
# 
# group4 <- long_soil %>% 
#   filter(Ion != "Ca", Ion != "K", Ion != "Mg", Ion != "NH4_N", 
#          Ion != "Mn", Ion != "Fe", Ion != "NO3_N",
#          Ion !="B", Ion !=  "Cd", Ion !=   "Cu",Ion != "Pb")
# group4
# 
# 
# boxplot1 <-  ggplot(group1, aes(x=Ion, y= Value, fill = treatment))+
#   geom_boxplot()
# boxplot1
# 
# boxplot2 <- ggplot(group2, aes(x=Ion, y= Value, fill = treatment))+
#   geom_boxplot()
# boxplot2
# 
# 
# boxplot3 <-  ggplot(group3, aes(x=Ion, y= Value, fill = treatment))+
#   geom_boxplot()
# boxplot3
# 
# boxplot4 <-  ggplot(group4, aes(x=Ion, y= Value, fill = treatment))+
#   geom_boxplot()
# boxplot4
# 
# ggarrange(boxplot1, boxplot2, boxplot3, boxplot4)

# PCA1: whole dataset ####
pca1_data <- soil %>% 
              select(NO3_N, NH4_N, Mg, K, P, Ca, Fe, Mn)

#Create dataframe that has all the sample information stored
var1 <- soil %>% 
  select(sample_n, site, type, treatment, pre_or_post, position) %>% 
  mutate(site = as.factor(site),
         type = as.factor(type),
         treatment = as.factor(treatment),
         pre_or_post = as.factor(pre_or_post),
         position = as.factor(position))

pca1 <- rda(pca1_data, scale = TRUE)
# summary(pca1, scaling = 1)

# PCA1 Broken stick model and Kaiser-Guttman criterion to select PC axes ####
# ev <- pca1$CA$eig #eigenvalues
# ev[ev>mean(ev)] #apply Kaiser-Guttman criterion to select axes: only the axes above the mean computed should be interpreted
# 
# #Broken stick model
# n <- length(ev)
# bsm <- data.frame(j=seq(1:n), p=0)
# for (i in 2:n) {
#   bsm$p[i] = bsm$p[i-1] + (1/(n+1-i))
# }
# 
# #Comparison of PCA result with broken stick model
# par(mfrow = c(2,1))
# barplot(ev, main = "Eigenvalues", col = "bisque", las=2)
# abline(h=mean(ev), col = "red")
# legend("topright", "Average eigenvalue", lwd = 1, col = 2, bty = "n")
# barplot(t(cbind(100*ev/sum(ev), bsm$p[n:1])), beside = TRUE, main = "% variance", col=c("bisque",2), las=2)
# 
# #plots using biplot.rda
# 
# # par(mfrow = c(1,2))
# # biplot.rda(pca1, scaling = 1, main = "PCA - Scaling 1")
# # biplot.rda(pca1, scaling = 2, main = "PCA - Scaling 2")
# #scaling1 = angles among descriptors meaningless, shows which samples have values clustered together
# #Scaling 2 distances among descriptors meaningless, shows which var1 are grouped together

# PCA2: only post-carcass drop ####

pca2_data <- soil %>% 
            filter(
              pre_or_post != "Baseline"
            )
rownames(pca2_data) <- pca2_data$sample_n
pca2_data <- pca2_data %>% mutate(sample_n = NULL, pre_or_post = NULL, site_type_treat = NULL)
pcapost <- pca2_data %>% 
           select(NO3_N, NH4_N, Ca, Mg, K, P, Fe, Mn)
var2 <- pca2_data %>% 
           select(site, treatment, type, position) %>% 
           mutate(site = as.factor(site), treatment = as.factor(treatment), type = as.factor(type), 
                  position = as.factor(position))
pca2 <- rda(pcapost, scale = TRUE)

# summary(pca2, scaling = 1)

# PCA2 Broken stick model and Kaiser-Guttman criterion to select PC axes ####
# ev2 <- pca2$CA$eig #eigenvalues
# ev2[ev2>mean(ev2)] #apply Kaiser-Guttman criterion to select axes: only the axes above the mean computed should be interpreted
# 
# #Broken stick model
# n <- length(ev2)
# bsm <- data.frame(j=seq(1:n), p=0)
# for (i in 2:n) {
#   bsm$p[i] = bsm$p[i-1] + (1/(n+1-i))
# }
# 
# par(mfrow = c(2,1))
# barplot(ev2, main = "Eigenvalues", col = "bisque", las=1)
# abline(h=mean(ev2), col = "red")
# legend("topright", "Average eigenvalue", lwd = 1, col = 2, bty = "n")
# barplot(t(cbind(100*ev/sum(ev2), bsm$p[n:1])), beside = TRUE, main = "% variance", col=c("bisque",2), las=2)
# 
# PCA3: only inside plot data, post carcass drop ####
pca3_data <- soil %>%
  filter(
    pre_or_post != "Baseline",
    position == "Inside"
  )
rownames(pca3_data) <- pca3_data$sample_n
pca3_data <- pca3_data %>% mutate(sample_n = NULL, pre_or_post = NULL, site_type_treat = NULL, position = NULL)
pcatreat <- pca3_data %>%
  select(NO3_N, NH4_N, Ca, Mg, K, P, Fe, Mn)
var3 <- pca3_data %>%
  select(site, treatment, type) %>%
  mutate(site = as.factor(site),
         treatment = as.factor(treatment),
         type = as.factor(type),
         )
# View(pca3_data)

pca3 <- rda(pcatreat, scale = TRUE)
# summary(pca3, scaling = 1)

# PCA3 Broken stick model and Kaiser-Guttman criterion to select PC axes ####
# ev3 <- pca3$CA$eig #eigenvalues
# ev3[ev3>mean(ev3)] #apply Kaiser-Guttman criterion to select axes: only the axes above the mean computed should be interpreted
# 
# #Broken stick model
# n <- length(ev3)
# bsm <- data.frame(j=seq(1:n), p=0)
# for (i in 2:n) {
#   bsm$p[i] = bsm$p[i-1] + (1/(n+1-i))
# }
# 
# par(mfrow = c(2,1))
# barplot(ev3, main = "Eigenvalues", col = "bisque", las=1)
# abline(h=mean(ev3), col = "red")
# legend("topright", "Average eigenvalue", lwd = 1, col = 2, bty = "n")
# barplot(t(cbind(100*ev/sum(ev3), bsm$p[n:1])), beside = TRUE, main = "% variance", col=c("bisque",2), las=2)
# 
# Get PCA scores and save them into excel spreadsheet ####
# results1 <- as.data.frame(pca1$CA$v[,c(1,2,3)]) %>% 
#   mutate(ion = rownames(results1),.before = PC1,
#          pca_n = vector(mode = "character", length = 8),
#          pca_n = "pca1")
# results2 <- as.data.frame(pca2$CA$v[,c(1,2,3)]) %>% 
#   mutate(ion = rownames(results2),.before = PC1,
#          pca_n = vector(mode = "character", length = 8),
#          pca_n = "pca2")
# results3 <- as.data.frame(pca3$CA$v[,c(1,2,3)]) %>% 
#   mutate(ion = rownames(results3),.before = PC1,
#          pca_n = vector(mode = "character", length = 8),
#          pca_n = "pca3")
# 
# results_soilPCAs <- rbind(results1, results2, results3) #rbind is to stack the dfs on top of each other
#                                                         #column names must be the same
# writexl::write_xlsx(results_soilPCAs, path = "pca1_results.xlsx")


# Plots PCA1 ####
percent_explained1 <- format(round((100 * pca1$CA$eig/sum(pca1$CA$eig)), digits = 1), nsmall = 1, trim = TRUE)

colvec <- c("darkblue", "orange", "grey68")
attach(var1)

# Plots PCA1 Axis 1, 2, 3, grouped by site####
png(filename = "figures/vegsoil/soil/pca1_bysite.png",
    width = 1415,
    height = 501,
    units = "px")

par(mfrow = c(1, 3))
ordiplot(pca1, choices = c(1,2), type = "n", 
         main = "Soil nutrients grouped by site, Axis 1 and 2", 
         xlab = glue("PCA Axis 1 ({percent_explained1[1]}%)"),
         ylab = glue("PCA Axis 2 ({percent_explained1[2]}%)"),
         cex.lab=1.5, cex.axis=1.5, cex.main=1.5, cex.sub=1.5)
    points(pca1, choices = c(1,2), col = colvec[site])
    text(pca1, choices = c(1,2),display = "species", 
         scaling = 2, col = "red", cex = 1.5)
    ordihull(pca1, choices = c(1,2), groups = site, col = colvec)
    legend("topright", legend = levels(site), bty = "n",
           col = colvec, pch = 21, pt.bg = colvec, cex = 1.5)
#PCA1 Axis 1, 3 grouped by site
ordiplot(pca1, choices = c(1,3), type = "n", 
         main = "Soil nutrients grouped by site, Axis 1 and 3", 
         xlab = glue("PCA Axis 1 ({percent_explained1[1]}%)"),
         ylab = glue("PCA Axis 3 ({percent_explained1[3]}%)"),
         cex.lab=1.5, cex.axis=1.5, cex.main=1.5, cex.sub=1.5)
    points(pca1, choices = c(1,3), col = colvec[site])
    text(pca1, choices = c(1,3), display = "species", 
         scaling = 2, col = "red", cex = 1.5)
    ordihull(pca1, choices = c(1,3), groups = site, col = colvec)
    legend("topright", legend = levels(site), bty = "n",
           col = colvec, pch = 21, pt.bg = colvec, cex = 1.5)
    
#PCA1 Axis 2, 3 grouped by site
ordiplot(pca1, choices = c(2,3), type = "n", 
             main = "Soil nutrients grouped by site, Axis 2 and 3", 
             xlab = glue("PCA Axis 2 ({percent_explained1[2]}%)"),
             ylab = glue("PCA Axis 3 ({percent_explained1[3]}%)"),
         cex.lab=1.5, cex.axis=1.5, cex.main=1.5, cex.sub=1.5)
    points(pca1, choices = c(2,3), col = colvec[site])
    text(pca1, choices = c(2,3),display = "species", 
         scaling = 2, col = "red", cex = 1.5)
    ordihull(pca1, choices = c(2,3), groups = site, col = colvec)
    legend("topright", legend = levels(site), bty = "n",
           col = colvec, pch = 21, pt.bg = colvec, cex = 1.5)
dev.off()

# Plots PCA1 Axis 1, 2, 3, grouped by type####    
png(filename = "figures/vegsoil/soil/pca1_bytype.png",
    width = 1415,
    height = 501,
    units = "px")

par(mfrow = c(1, 3))

ordiplot(pca1, choices = c(1,2), type = "n", 
             main = "Soil nutrients grouped by type, Axis 1 and 2", 
             xlab = glue("PCA Axis 1 ({percent_explained1[1]}%)"),
             ylab = glue("PCA Axis 2 ({percent_explained1[2]}%)"),
         cex.lab=1.5, cex.axis=1.5, cex.main=1.5, cex.sub=1.5)
    points(pca1, choices = c(1,2), col = colvec[type])
    text(pca1, choices = c(1,2),display = "species", 
         scaling = 2, col = "red", cex = 1.5)
    ordihull(pca1, choices = c(1,2), groups = type, col = colvec)
    legend("topright", legend = levels(type), bty = "n",
           col = colvec, pch = 21, pt.bg = colvec, cex = 1.5)
#PCA1 type Axis 1, 3 grouped by type
ordiplot(pca1, choices = c(1,3), type = "n", 
             main = "Soil nutrients grouped by type, Axis 1 and 3", 
             xlab = glue("PCA Axis 1 ({percent_explained1[1]}%)"),
             ylab = glue("PCA Axis 3 ({percent_explained1[3]}%)"),
         cex.lab=1.5, cex.axis=1.5, cex.main=1.5, cex.sub=1.5)
    points(pca1, choices = c(1,3), col = colvec[type])
    text(pca1, choices = c(1,3),display = "species", 
         scaling = 2, col = "red", cex = 1.5)
    ordihull(pca1, choices = c(1,3), groups = type, col = colvec)
    legend("topright", legend = levels(type), bty = "n",
           col = colvec, pch = 21, pt.bg = colvec, cex = 1.5)
    
#PCA1 type Axis 2, 3 grouped by type
ordiplot(pca1, choices = c(2,3), type = "n", 
             main = "Soil nutrients grouped by type, Axis 2 and 3", 
             xlab = glue("PCA Axis 2 ({percent_explained1[2]}%)"),
             ylab = glue("PCA Axis 3 ({percent_explained1[3]}%)"),
         cex.lab=1.5, cex.axis=1.5, cex.main=1.5, cex.sub=1.5)
    points(pca1, choices = c(2,3), col = colvec[type])
    text(pca1, choices = c(2,3),display = "species", 
         scaling = 2, col = "red", cex = 1.5)
    ordihull(pca1, choices = c(2,3), groups = type, col = colvec)
    legend("topright", legend = levels(type), bty = "n",
           col = colvec, pch = 21, pt.bg = colvec, cex = 1.5)

dev.off()
# Plots PCA1 Axis 1, 2, 3, grouped by treatment####
png(filename = "figures/vegsoil/soil/pca1_bytreat.png",
    width = 1415,
    height = 501,
    units = "px")

par(mfrow = c(1, 3))
ordiplot(pca1, choices = c(1,2), type = "n", 
             main = "Soil nutrients grouped by treatment, Axis 1 and 2", 
             xlab = glue("PCA Axis 1 ({percent_explained1[1]}%)"),
             ylab = glue("PCA Axis 2 ({percent_explained1[2]}%)"),
         cex.lab=1.5, cex.axis=1.5, cex.main=1.5, cex.sub=1.5)
    points(pca1, choices = c(1,2), col = colvec[treatment])
    text(pca1, choices = c(1,2),display = "species", 
         scaling = 2, col = "red", cex = 1.5)
    ordihull(pca1, choices = c(1,2), groups = treatment, col = colvec)
    legend("topright", legend = levels(treatment), bty = "n",
           col = colvec, pch = 21, pt.bg = colvec, cex = 1.5)
#PCA1 type Axis 1, 3 grouped by treatment
ordiplot(pca1, choices = c(1,3), type = "n", 
             main = "Soil nutrients grouped by treatment, Axis 1 and 3", 
             xlab = glue("PCA Axis 1 ({percent_explained1[1]}%)"),
             ylab = glue("PCA Axis 3 ({percent_explained1[3]}%)"),
         cex.lab=1.5, cex.axis=1.5, cex.main=1.5, cex.sub=1.5)
    points(pca1, choices = c(1,3), col = colvec[treatment])
    text(pca1, choices = c(1,3),display = "species", 
         scaling = 2, col = "red", cex = 1.5)
    ordihull(pca1, choices = c(1,3), groups = treatment, col = colvec)
    legend("topright", legend = levels(treatment), bty = "n",
           col = colvec, pch = 21, pt.bg = colvec, cex = 1.5)
    
#PCA1 type Axis 2, 3 grouped by treatment
ordiplot(pca1, choices = c(2,3), type = "n", 
             main = "Soil nutrients grouped by treatment, Axis 2 and 3", 
             xlab = glue("PCA Axis 2 ({percent_explained1[2]}%)"),
             ylab = glue("PCA Axis 3 ({percent_explained1[3]}%)"),
         cex.lab=1.5, cex.axis=1.5, cex.main=1.5, cex.sub=1.5)
    points(pca1, choices = c(2,3), col = colvec[treatment])
    text(pca1, choices = c(2,3),display = "species", 
         scaling = 2, col = "red", cex = 1.5)
    ordihull(pca1, choices = c(2,3), groups = treatment, col = colvec)
    legend("topright", legend = levels(treatment), bty = "n",
           col = colvec, pch = 21, pt.bg = colvec, cex = 1.5)
dev.off()
# Plots PCA1 Axis 1, 2, 3, grouped by position####
png(filename = "figures/vegsoil/soil/pca1_byposition.png",
    width = 1415,
    height = 501,
    units = "px")

par(mfrow = c(1, 3))
ordiplot(pca1, choices = c(1,2), type = "n", 
             main = "Soil nutrients inside vs outside plots, Axis 1 and 2", 
             xlab = glue("PCA Axis 1 ({percent_explained1[1]}%)"),
             ylab = glue("PCA Axis 2 ({percent_explained1[2]}%)"),
         cex.lab=1.5, cex.axis=1.5, cex.main=1.5, cex.sub=1.5)
    points(pca1, choices = c(1,2), col = colvec[position])
    text(pca1, choices = c(1,2),display = "species", 
         scaling = 2, col = "red", cex = 1.5)
    ordihull(pca1, choices = c(1,2), groups = position, col = colvec)
    legend("topright", legend = levels(position), bty = "n",
           col = colvec, pch = 21, pt.bg = colvec, cex = 1.5)
#PCA1 type Axis 1, 3 grouped by position
ordiplot(pca1, choices = c(1,3), type = "n", 
             main = "Soil nutrients inside vs outside plots, Axis 1 and 3", 
             xlab = glue("PCA Axis 1 ({percent_explained1[1]}%)"),
             ylab = glue("PCA Axis 3 ({percent_explained1[3]}%)"),
         cex.lab=1.5, cex.axis=1.5, cex.main=1.5, cex.sub=1.5)
    points(pca1, choices = c(1,3), col = colvec[position])
    text(pca1, choices = c(1,3),display = "species", 
         scaling = 2, col = "red", cex = 1.5)
    ordihull(pca1, choices = c(1,3), groups = position, col = colvec)
    legend("topright", legend = levels(position), bty = "n",
           col = colvec, pch = 21, pt.bg = colvec, cex = 1.5)
    
#PCA1 type Axis 2, 3 grouped by position
ordiplot(pca1, choices = c(2,3), type = "n", 
             main = "Soil nutrients inside vs outside plots, Axis 2 and 3", 
             xlab = glue("PCA Axis 2 ({percent_explained1[2]}%)"),
             ylab = glue("PCA Axis 3 ({percent_explained1[3]}%)"),
         cex.lab=1.5, cex.axis=1.5, cex.main=1.5, cex.sub=1.5)
    points(pca1, choices = c(2,3), col = colvec[position])
    text(pca1, choices = c(2,3),display = "species", 
         scaling = 2, col = "red", cex = 1.5)
    ordihull(pca1, choices = c(2,3), groups = position, col = colvec)
    legend("topright", legend = levels(position), bty = "n",
           col = colvec, pch = 21, pt.bg = colvec, cex = 1.5)
dev.off()
# Plots PCA1 Axis 1, 2, 3, grouped by pre_or_post####
png(filename = "figures/vegsoil/soil/pca1_preorpost.png",
    width = 1415,
    height = 501,
    units = "px")
par(mfrow = c(1, 3))
ordiplot(pca1, choices = c(1,2), type = "n", 
             main = "Soil nutrients: baseline vs post data, Axis 1 and 2", 
             xlab = glue("PCA Axis 1 ({percent_explained1[1]}%)"),
             ylab = glue("PCA Axis 2 ({percent_explained1[2]}%)"),
         cex.lab=1.5, cex.axis=1.5, cex.main=1.5, cex.sub=1.5)
    points(pca1, choices = c(1,2), col = colvec[pre_or_post])
    text(pca1, choices = c(1,2),display = "species", 
         scaling = 2, col = "red", cex = 1.5)
    ordihull(pca1, choices = c(1,2), groups = pre_or_post, col = colvec)
    legend("topright", legend = levels(pre_or_post), bty = "n",
           col = colvec, pch = 21, pt.bg = colvec, cex = 1.5)
#PCA1 type Axis 1, 3 grouped by pre_or_post
ordiplot(pca1, choices = c(1,3), type = "n", 
             main = "Soil nutrients: baseline vs post data, Axis 1 and 3", 
             xlab = glue("PCA Axis 1 ({percent_explained1[1]}%)"),
             ylab = glue("PCA Axis 3 ({percent_explained1[3]}%)"),
         cex.lab=1.5, cex.axis=1.5, cex.main=1.5, cex.sub=1.5)
    points(pca1, choices = c(1,3), col = colvec[pre_or_post])
    text(pca1, choices = c(1,3),display = "species", 
         scaling = 2, col = "red", cex = 1.5)
    ordihull(pca1, choices = c(1,3), groups = pre_or_post, col = colvec)
    legend("topright", legend = levels(pre_or_post), bty = "n",
           col = colvec, pch = 21, pt.bg = colvec, cex = 1.5)
    
#PCA1 type Axis 2, 3 grouped by pre_or_post
ordiplot(pca1, choices = c(2,3), type = "n", 
             main = "Soil nutrients: baseline vs post data, Axis 2 and 3", 
             xlab = glue("PCA Axis 2 ({percent_explained1[2]}%)"),
             ylab = glue("PCA Axis 3 ({percent_explained1[3]}%)"),
         cex.lab=1.5, cex.axis=1.5, cex.main=1.5, cex.sub=1.5)
    points(pca1, choices = c(2,3), col = colvec[pre_or_post])
    text(pca1, choices = c(2,3),display = "species", 
         scaling = 2, col = "red", cex = 1.5)
    ordihull(pca1, choices = c(2,3), groups = pre_or_post, col = colvec)
    legend("topright", legend = levels(pre_or_post), bty = "n",
           col = colvec, pch = 21, pt.bg = colvec, cex = 1.5)
dev.off()  
# Plots for PCA2 ####
percent_explained2 <- format(round((100 * pca2$CA$eig/sum(pca2$CA$eig)), digits = 1), nsmall = 1, trim = TRUE)
colvec <- c("darkblue", "orange", "grey68")
    
detach(var1)
attach(var2)
# Plots pca2 Axis 1, 2, 3, grouped by site####
png(filename = "figures/vegsoil/soil/pca2_bysite.png",
    width = 1415,
    height = 501,
    units = "px")

par(mfrow = c(1, 3))
ordiplot(pca2, choices = c(1,2), type = "n", 
             main = "Soil nutrients grouped by site, Axis 1 and 2", 
             xlab = glue("PCA Axis 1 ({percent_explained2[1]}%)"),
             ylab = glue("PCA Axis 2 ({percent_explained2[2]}%)"),
         cex.lab=1.5, cex.axis=1.5, cex.main=1.5, cex.sub=1.5)
    points(pca2, choices = c(1,2), col = colvec[site])
    text(pca2, choices = c(1,2),display = "species", 
         scaling = 2, col = "red", cex = 1.5)
    ordihull(pca2, choices = c(1,2), groups = site, col = colvec)
    legend("topright", legend = levels(site), bty = "n",
           col = colvec, pch = 21, pt.bg = colvec, cex = 1.5)
# Plots pca2 Axis 1, 3 grouped by site
ordiplot(pca2, choices = c(1,3), type = "n", 
             main = "Soil nutrients grouped by site, Axis 1 and 3", 
             xlab = glue("PCA Axis 1 ({percent_explained2[1]}%)"),
             ylab = glue("PCA Axis 3 ({percent_explained2[3]}%)"),
         cex.lab=1.5, cex.axis=1.5, cex.main=1.5, cex.sub=1.5)
    points(pca2, choices = c(1,3), col = colvec[site])
    text(pca2, choices = c(1,3),display = "species", 
         scaling = 2, col = "red", cex = 1.5)
    ordihull(pca2, choices = c(1,3), groups = site, col = colvec)
    legend("topright", legend = levels(site), bty = "n",
           col = colvec, pch = 21, pt.bg = colvec, cex = 1.5)
    
# Plots pca2 Axis 2, 3 grouped by site
ordiplot(pca2, choices = c(2,3), type = "n", 
             main = "Soil nutrients grouped by site, Axis 2 and 3", 
             xlab = glue("PCA Axis 2 ({percent_explained2[2]}%)"),
             ylab = glue("PCA Axis 3 ({percent_explained2[3]}%)"),
         cex.lab=1.5, cex.axis=1.5, cex.main=1.5, cex.sub=1.5)
    points(pca2, choices = c(2,3), col = colvec[site])
    text(pca2, choices = c(2,3),display = "species", 
         scaling = 2, col = "red", cex = 1.5)
    ordihull(pca2, choices = c(2,3), groups = site, col = colvec)
    legend("topright", legend = levels(site), bty = "n",
           col = colvec, pch = 21, pt.bg = colvec, cex = 1.5)
dev.off()

# Plots pca2 Axis 1, 2, 3, grouped by type####   
png(filename = "figures/vegsoil/soil/pca2_bytype.png",
    width = 1415,
    height = 501,
    units = "px")
par(mfrow = c(1, 3))
ordiplot(pca2, choices = c(1,2), type = "n", 
             main = "Soil nutrients grouped by type, Axis 1 and 2", 
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
             main = "Soil nutrients grouped by type, Axis 1 and 3", 
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
             main = "Soil nutrients grouped by type, Axis 2 and 3", 
             xlab = glue("PCA Axis 2 ({percent_explained2[2]}%)"),
             ylab = glue("PCA Axis 3 ({percent_explained2[3]}%)"),
         cex.lab=1.5, cex.axis=1.5, cex.main=1.5, cex.sub=1.5)
    points(pca2, choices = c(2,3), col = colvec[type])
    text(pca2, choices = c(2,3),display = "species", 
         scaling = 2, col = "red", cex = 1.5)
    ordihull(pca2, choices = c(2,3), groups = type, col = colvec)
    legend("topleft", legend = levels(type), bty = "n",
           col = colvec, pch = 21, pt.bg = colvec, cex = 1.5)
dev.off()    
# Plots pca2 Axis 1, 2, 3, grouped by treatment####
png(filename = "figures/vegsoil/soil/pca2_bytreat.png",
    width = 1415,
    height = 501,
    units = "px")

par(mfrow = c(1, 3))
ordiplot(pca2, choices = c(1,2), type = "n", 
             main = "Soil nutrients grouped by treatment, Axis 1 and 2", 
             xlab = glue("PCA Axis 1 ({percent_explained2[1]}%)"),
             ylab = glue("PCA Axis 2 ({percent_explained2[2]}%)"),
         cex.lab=1.5, cex.axis=1.5, cex.main=1.5, cex.sub=1.5)
    points(pca2, choices = c(1,2), col = colvec[treatment])
    text(pca2, choices = c(1,2),display = "species", 
         scaling = 2, col = "red", cex = 1.5)
    ordihull(pca2, choices = c(1,2), groups = treatment, col = colvec)
    legend("topleft", legend = levels(treatment), bty = "n",
           col = colvec, pch = 21, pt.bg = colvec, cex = 1.5)
# Plots pca2 Axis 1, 3 grouped by treatment
ordiplot(pca2, choices = c(1,3), type = "n", 
             main = "Soil nutrients grouped by treatment, Axis 1 and 3", 
             xlab = glue("PCA Axis 1 ({percent_explained2[1]}%)"),
             ylab = glue("PCA Axis 3 ({percent_explained2[3]}%)"),
         cex.lab=1.5, cex.axis=1.5, cex.main=1.5, cex.sub=1.5)
    points(pca2, choices = c(1,3), col = colvec[treatment])
    text(pca2, choices = c(1,3),display = "species", 
         scaling = 2, col = "red", cex = 1.5)
    ordihull(pca2, choices = c(1,3), groups = treatment, col = colvec)
    legend("topleft", legend = levels(treatment), bty = "n",
           col = colvec, pch = 21, pt.bg = colvec, cex = 1.5)
    
# Plots pca2 Axis 2, 3 grouped by treatment
ordiplot(pca2, choices = c(2,3), type = "n", 
             main = "Soil nutrients grouped by treatment, Axis 2 and 3", 
             xlab = glue("PCA Axis 2 ({percent_explained2[2]}%)"),
             ylab = glue("PCA Axis 3 ({percent_explained2[3]}%)"),
         cex.lab=1.5, cex.axis=1.5, cex.main=1.5, cex.sub=1.5)
    points(pca2, choices = c(2,3), col = colvec[treatment])
    text(pca2, choices = c(2,3),display = "species", 
         scaling = 2, col = "red", cex = 1.5)
    ordihull(pca2, choices = c(2,3), groups = treatment, col = colvec)
    legend("topleft", legend = levels(treatment), bty = "n",
           col = colvec, pch = 21, pt.bg = colvec, cex = 1.5)
    dev.off()
# Plots pca2 Axis 1, 2, 3, grouped by position####
png(filename = "figures/vegsoil/soil/pca2_byposition.png",
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
             main = "Soil nutrients inside vs outside plots, Axis 2 and 3", 
             xlab = glue("PCA Axis 2 ({percent_explained2[2]}%)"),
             ylab = glue("PCA Axis 3 ({percent_explained2[3]}%)"),
             cex.lab=1.5, cex.axis=1.5, cex.main=1.5, cex.sub=1.5)
    points(pca2, choices = c(2,3), col = colvec[position])
    text(pca2, choices = c(2,3),display = "species", 
         scaling = 2, col = "red", cex = 1.5)
    ordihull(pca2, choices = c(2,3), groups = position, col = colvec)
    legend("topright", legend = levels(position), bty = "n",
           col = colvec, pch = 21, pt.bg = colvec, cex = 1.5)
dev.off()
# Plots for PCA3  ####    
percent_explained3<- format(round((100 * pca3$CA$eig/sum(pca3$CA$eig)), digits = 1), nsmall = 1, trim = TRUE)
colvec <- c("darkblue", "orange", "grey68")
    
detach(var2)
attach(var3)
# Plots pca3 Axis 1, 2, 3, grouped by site####
png(filename = "figures/vegsoil/soil/pca3_bysite.png",
    width = 1415,
    height = 501,
    units = "px")

par(mfrow = c(1, 3))
ordiplot(pca3, choices = c(1,2), type = "n", 
             main = "Soil nutrients grouped by site, Axis 1 and 2", 
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
             main = "Soil nutrients grouped by site, Axis 1 and 3", 
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
             main = "Soil nutrients grouped by site, Axis 2 and 3", 
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
# Plots pca3 Axis 1, 2, 3, grouped by type####  
png(filename = "figures/vegsoil/soil/pca3_bytype.png",
    width = 1415,
    height = 501,
    units = "px")

par(mfrow = c(1, 3))

ordiplot(pca3, choices = c(1,2), type = "n", 
             main = "Soil nutrients grouped by type, Axis 1 and 2", 
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
             main = "Soil nutrients grouped by type, Axis 1 and 3", 
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
             main = "Soil nutrients grouped by type, Axis 2 and 3", 
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
# Plots pca3 Axis 1, 2, 3, grouped by treatment####
png(filename = "figures/vegsoil/soil/pca3_bytreat.png",
    width = 1415,
    height = 501,
    units = "px")

par(mfrow = c(1, 3))

ordiplot(pca3, choices = c(1,2), type = "n", 
             main = "Soil nutrients grouped by treatment, Axis 1 and 2", 
             xlab = glue("PCA Axis 1 ({percent_explained3[1]}%)"),
             ylab = glue("PCA Axis 2 ({percent_explained3[2]}%)"),
         cex.lab=1.5, cex.axis=1.5, cex.main=1.5, cex.sub=1.5)
    points(pca3, choices = c(1,2), col = colvec[treatment])
    text(pca3, choices = c(1,2),display = "species", 
         scaling = 2, col = "red", cex= 1.5)
    ordihull(pca3, choices = c(1,2), groups = treatment, col = colvec)
    legend("topright", legend = levels(treatment), bty = "n",
           col = colvec, pch = 21, pt.bg = colvec, cex= 1.5)
# Plots pca3 Axis 1, 3 grouped by treatment
ordiplot(pca3, choices = c(1,3), type = "n", 
             main = "Soil nutrients grouped by treatment, Axis 1 and 3", 
             xlab = glue("PCA Axis 1 ({percent_explained3[1]}%)"),
             ylab = glue("PCA Axis 3 ({percent_explained3[3]}%)"),
         cex.lab=1.5, cex.axis=1.5, cex.main=1.5, cex.sub=1.5)
    points(pca3, choices = c(1,3), col = colvec[treatment])
    text(pca3, choices = c(1,3),display = "species", 
         scaling = 2, col = "red", cex= 1.5)
    ordihull(pca3, choices = c(1,3), groups = treatment, col = colvec)
    legend("topright", legend = levels(treatment), bty = "n",
           col = colvec, pch = 21, pt.bg = colvec, cex= 1.5)
    
# Plots pca3 Axis 2, 3 grouped by treatment
ordiplot(pca3, choices = c(2,3), type = "n", 
             main = "Soil nutrients grouped by treatment, Axis 2 and 3", 
             xlab = glue("PCA Axis 2 ({percent_explained3[2]}%)"),
             ylab = glue("PCA Axis 3 ({percent_explained3[3]}%)"),
         cex.lab=1.5, cex.axis=1.5, cex.main=1.5, cex.sub=1.5)
    points(pca3, choices = c(2,3), col = colvec[treatment])
    text(pca3, choices = c(2,3),display = "species", 
         scaling = 2, col = "red", cex= 1.5)
    ordihull(pca3, choices = c(2,3), groups = treatment, col = colvec)
    legend("topright", legend = levels(treatment), bty = "n",
           col = colvec, pch = 21, pt.bg = colvec, cex= 1.5)
dev.off()