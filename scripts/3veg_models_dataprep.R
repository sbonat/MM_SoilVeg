#packages ####
# General data reading and manipulation
library(readxl)
library(tidyverse)
library(lubridate)
library(broom) #useful for tidying model outputs
library(broom.mixed) #for mixed models
library(broom.helpers) #labelling etc
library(labelled)

# Graphics
library(ggplot2)
library(ggpubr)
library(sjPlot)
library(lares)

# Statistics
library(nlme)
library(lme4)
library(MuMIn)
library(itsadug) #has acf_resid() formula

#Model validation and other
library(papeR)
library(performance)
library(predictmeans)

# Load data ####
veg0 <- read_excel("data/veg_clean.xlsx")
veg <- veg0 %>% 
  mutate(pre_or_post = gsub("Pre", "Baseline", pre_or_post),
         plot_id = paste(site,type,treatment, sep = ""),
         plot_id = gsub("Nil", "NIL", plot_id),
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
         date = ymd(date)
         ) %>% 
  filter(pre_or_post != "Year1") %>% 
  rename(debris = deadmaterial)%>%
  mutate_at(c(13:19), ~.x/100) #Divide response columns by 100 to get proportions

# Divide post by pre data
veg_pre <- veg %>% 
  filter(pre_or_post == "Baseline")  %>% 
  arrange(plot_id, position, transect, quadrat)    
# Select the variables which will need to be merged back into the dataframe later
veg_var <- veg_pre %>% 
  select(c(1:12))
# Then select the columns to divide
veg_pre <- veg_pre %>%
  select(c(13:19))

# Select data from Post carcass drop surveys and arrange by plot_id and position
veg_post <- veg %>%
  filter(pre_or_post == "Post")  %>% 
  arrange(plot_id, position, transect, quadrat) %>%
  select(c(13:19))

#Subtract Post - Pre. In this case it's better to do this, as we have proportion data so dividing
#like was done for the soil data would not be ideal

veg_subtracted <- veg_post - veg_pre

# Merge back with variables
veg1 <- cbind(veg_var, veg_subtracted)

veg1 <-  veg1 %>% 
        mutate(site = as.factor(site),
         type = as.factor(type),
         treatment = as.factor(treatment),
         pre_or_post = NULL,
         date = NULL,
         transect = as.factor(transect),
         quadrat = as.factor(quadrat),
         position = as.factor(position))

# Tidy up environment
rm(veg_subtracted,veg_pre, veg_post, veg_var, veg, veg0)

# Data explovaluen plots ####
veg_nested <- veg1 %>% 
  pivot_longer(cols = c(11:17),
               names_to = "response",
               values_to = "value") %>% 
  group_by(response) %>% 
  nest()

# Histograms
par(mfrow = c(3, 3))
veg_nested1 <- veg_nested %>% 
  mutate(hist = map(data, ~ hist(.x$value, 
                                 xlab = "Change in cover (Post-Baseline)",
                                 ylab = "Frequency",
                                 main = paste(response)))) %>% 
  select(-hist)

# Cleveland dotplots
par(mfrow = c(3, 3))
veg_nested1 <- veg_nested %>% 
  mutate(dotplots = map(data, ~dotchart(.x$value, 
                                        xlab = "Change in cover (Post-Baseline)",
                                        ylab = "Order of observations",
                                        main = paste(response)))) %>% 
  select(-dotplots)

# Plot each variable against explanatory variables. First, response vs veg height
par(mfrow = c(3, 3))
veg_nested <- veg_nested %>%
  mutate(plots_height = map(data, ~plot(x = .x$value, y = .x$height,
                                        xlab = "Change (Post-Baseline)",
                                        ylab = "Height (cm)",
                                        main = paste(response)))) %>% 
  select(-plots_height)

# Second, response vs moisture
par(mfrow = c(3, 3))
veg_nested <- veg_nested %>%
  mutate(plots_moisture = map(data, ~plot(x = .x$value, y = .x$moisture,
                                          xlab = "Change (Post-Baseline)",
                                          ylab = "Moisture (%v)",
                                          main = paste(response)))) %>%
  select(-plots_moisture)

# Reset plot pane
dev.off()

#At this stage we can select out mossfern and shrub, as they have very little data
rm(veg_nested1)
veg_nested <- veg_nested %>% 
  filter(response != "mossfern",
         response != "shrub")

veg1 <- veg1 %>% 
  select(-mossfern, -shrub)

#Double check correlations among variables, but do it for each response separately

#Debris
corr_data1 <- veg1 %>% 
  dplyr::select(site, type, treatment, position, transect, quadrat, 
                height, moisture, debris)
corr_cross(corr_data1, max_pvalue = 0.05, top = 10)

#Bareground
corr_data2 <- veg1 %>% 
  dplyr::select(site, type, treatment, position, transect, quadrat, 
                height, moisture, bareground)
corr_cross(corr_data2, max_pvalue = 0.05, top = 10)

#Cover
corr_data3 <- veg1 %>% 
  dplyr::select(site, type, treatment, position, transect, quadrat, 
                height, moisture, bareground)
corr_cross(corr_data3, max_pvalue = 0.05, top = 10)

#Forb
corr_data4 <- veg1 %>% 
  dplyr::select(site, type, treatment, position, transect, quadrat, 
                height, moisture, bareground)
corr_cross(corr_data4, max_pvalue = 0.05, top = 10)

#Grasses
corr_data5 <- veg1 %>% 
  dplyr::select(site, type, treatment, position, transect, quadrat, 
                height, moisture, bareground)
corr_cross(corr_data5, max_pvalue = 0.05, top = 10)

#Site might cause some collinearity issues.
#Not interested in its effect, should really be a random effect
#Site, transect, quadrat to go as random effects
#Moisture and height are highly correlated with site, so should be taken into account by random effect
#Position might cause collinearity

rm(corr_data1, corr_data2, corr_data3, corr_data4, corr_data5)

#Further data manipulation: inside vs outside plots cause singularity issues####
#This bit of code might need to be tidied

veg1 <- veg1 %>% 
  mutate(pre_or_post = NULL,
         date = NULL)%>% 
  filter(position == "Inside") %>% #position is causing singularity issues so just focus on Inside plot
  select(-position)

plot_coord <- read_excel("data/site_coordinates.xlsx")

veg1 <- merge(plot_coord, veg1, by = "plot_id")

rm(plot_coord)

veg1 <- veg1 %>% 
  select(-bareground, -debris)

#Get soil probe data for N and P which will be useful for veg data analysis
soil <- read_excel("data/soil_analysis.xlsx")

soil_v <- soil %>% 
  filter(position == "Inside") %>% 
  select(plot_id, sample_n, site, treatment, type, inorg_N, P, K, Mg) #select the nutrients most relevant for plant growth

veg2 <- merge(soil_v, veg1, by = c("plot_id", "site", "treatment", "type")) 

#writexl::write_xlsx(veg2, "data/veg_analysis.xlsx")
rm(soil_v, veg2)