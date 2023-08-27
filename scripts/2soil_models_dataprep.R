# Packages ####

# General data reading and manipulation
library(readxl)
library(tidyverse)
library(lubridate)

# Graphics
library(ggplot2)
library(ggpubr)
library(lares) #cross correlations

# Load data####
soil0 <- read_excel("data/soil_probes_clean.xlsx")
soil <- soil0 %>% 
  rename(timing = pre_or_post) %>% 
  mutate(treatment = gsub("HE", "Open", treatment),
         treatment = gsub("Nil", "Open", treatment),
         treatment = gsub("VE", "Vertebrate exclusion", treatment),
         treatment = gsub("IE", "Invertebrate exclusion", treatment),
         type = gsub("MM", "Mass mortality", type),
         type = gsub("SC", "Single carcass", type),
         type = gsub("C", "Control", type),
         site = gsub("CR", "Site 1", site),
         site = gsub("MS", "Site 3", site),
         site = gsub("WP", "Site 2", site),
         timing = gsub("Pre", "Baseline", timing),
         treatment = gsub("NIL", "Open", treatment),
         B = NULL, Cd = NULL, Cu = NULL, Pb = NULL, Zn = NULL, Al = NULL, S = NULL,
         '#anion' = NULL, '#cation' = NULL, retrieval_date = NULL,
         treatment = as.factor(treatment),
         site = as.factor(site),
         type = as.factor(type),
         position = as.factor(position),
         site_type_treat = as.factor(site_type_treat)
  ) %>% 
  filter(timing != "Yr1", # selecting out Yr1 results as not relevant
         sample_n != "5", # sample #5 and 6 seem to have been entered wrong so have been removed from the dataset
         sample_n != "6") %>%  
  rename("plot_id" = "site_type_treat")
soil <- soil %>% 
  mutate(
    timing = as.factor(timing))%>% 
  mutate(
    inorg_N = NH4_N + NO3_N, .after = burial_date, #Added together NO3 and NH4 to look at them together as Inorganic Nitrogen
    NH4_N = NULL, NO3_N = NULL)

#Add coordinates to plots
plot_coord <- read_excel("data/site_coordinates.xlsx")
soil <-  merge(plot_coord, soil, by = "plot_id")
rm(plot_coord)

# Data manipulation####
# Dividing Post/Pre data: this should normalize data, and it should take into account variation between sites.
# It will also mean I have less variables to include in models and the time variable
# has to be interpreted when looking at estimates as a proportional increase/decrease

# First, select data from Baseline surveys and arrange by plot_id and position
soil_pre <- soil %>% 
  filter(timing == "Baseline")  %>% 
  arrange(plot_id, position)    
#Select the variables which will need to be merged back into the dataframe later
soil_var <- soil_pre %>% 
  select(c(1:8))
#Then select the columns to divide
soil_pre <- soil_pre %>%
  select(inorg_N, Ca, Mg, Mn, Fe, K, P)

#Select data from Post carcass drop surveys and arrange by plot_id and position
soil_post <- soil %>%
  filter(timing == "Post")  %>% 
  arrange(plot_id, position) %>%
  select(inorg_N, Ca, Mg, Mn, Fe, K, P)

soil_divided <- soil_post/soil_pre

# Merge back with variables
soil1 <- cbind(soil_var, soil_divided)

#writexl::write_xlsx(soil1, "data/soil_analysis.xlsx")
#Tidy up environment
rm(soil_divided, soil_pre, soil_post, soil_var)

#Nest dataframe and run histograms and dotcharts
soil_nested <- soil1 %>% 
  pivot_longer(cols = c(9:15),
               names_to = "nutrient",
               values_to = "ratio"
  ) %>% 
  group_by(nutrient) %>% 
  nest() 

# Data exploration plots####
par(mfrow = c(3,3))

soil_nested <-  soil_nested%>% 
  mutate(hist1 = map(data, ~hist(.x$ratio, 
                                 xlab = "Nutrient concentration", 
                                 ylab = "Frequency", main = paste(nutrient)))) %>% 
  select(-hist1)


par(mfrow = c(3,3))
soil_nested <-  soil_nested%>% 
  mutate(dotplots = map(data, ~dotchart(.x$ratio,
                                        xlab = "Nutrient concentration", 
                                        ylab = "Order of observations", main = paste(nutrient)))) %>% 
  select(-dotplots)

dev.off()

# can expect some homogeneity of variance issues

# Check correlations among variables

#Inorg_N
corr_data1 <- soil1 %>% 
  dplyr::select(site, type, treatment, position, 
                inorg_N)

corr_cross(corr_data1, max_pvalue = 0.05, top = 20)

#Ca
corr_data2 <- soil1 %>% 
  dplyr::select(site, type, treatment, position,  
                Ca)
corr_cross(corr_data2, max_pvalue = 0.05, top = 20)

#Mg
corr_data3 <- soil1 %>% 
  dplyr::select(site, type, treatment, position, 
                Mg)
corr_cross(corr_data3, max_pvalue = 0.05, top = 20)

#Fe
corr_data4 <- soil1 %>% 
  dplyr::select(site, type, treatment, position,  
                Fe)
corr_cross(corr_data4, max_pvalue = 0.05, top = 20)

#Mn
corr_data5 <- soil1 %>% 
  dplyr::select(site, type, treatment, position,  
                Mn)
corr_cross(corr_data5, max_pvalue = 0.05, top = 20)

#K
corr_data6 <- soil1 %>% 
  dplyr::select(site, type, treatment, position,  
                K)
corr_cross(corr_data6, max_pvalue = 0.05, top = 20)

#P
corr_data7 <- soil1 %>% 
  dplyr::select(site, type, treatment, position,  
                P)
corr_cross(corr_data7, max_pvalue = 0.05, top = 20)