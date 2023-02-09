# Packages ####

# General data reading and manipulation
library(readxl)
library(tidyverse)
library(lubridate)
library(broom) #useful for tidying model outputs
library(broom.mixed) #for mixed models
library(broom.helpers) #labelling etc

# Graphics
library(ggplot2)
library(ggpubr)
library(lares)

# Statistics
library(nlme)
library(lme4)
library(MuMIn)

#Model validation and other
library(predictmeans)
library(emmeans) #calculate marginal means and plot


# Load data####
soil <- read_excel("data/soil_analysis.xlsx")

soil1 <- soil %>% 
  select(plot_id, sample_n, site, type, treatment, position, Fe) %>% 
  mutate(site = as.factor(site),
         type = as.factor(type),
         treatment = as.factor(treatment),
         position = as.factor(position)
  )
#See the data prep script within this project for all the data exploration

#iron ####

# Formulas. f1a, f1b, f1 are the formulas including all terms possible. 
# These are used to find the optimal variance structure, then to find the optimal random structure
# The rest of the formulas will be updated later when the above has been accomplished

f1a <- Fe ~ treatment + type + position + treatment*type + type*position + treatment*position

#Log transform data
f1b <- log(Fe) ~ treatment + type + position + treatment*type + type*position + treatment*position

# First, run a gls, then check residuals and ACF
m1 <-  gls(f1a,
           method = "REML",
           data = soil1)

residplot(m1, newwd=F)

#problems with homogeneity of variance and autocorrelation

# Now, see if adding a weight to the model improves its fit. Try with different variables

m1.v1 <-  gls(f1a,
              method = "REML",
              weights = varIdent(form = ~1|site),
              data = soil1)
m1.v2 <-  gls(f1a,
              method = "REML",
              weights = varIdent(form = ~1|type),
              data = soil1)
m1.v3 <-  gls(f1a,
              method = "REML",
              weights = varIdent(form = ~1|treatment),
              data = soil1)
m1.v4 <-  gls(f1a,
              method = "REML",
              weights = varIdent(form = ~1|position),
              data = soil1)

MuMIn::AICc(m1, m1.v1, m1.v2, m1.v3, m1.v4)
rm(m1, m1.v2, m1.v3, m1.v4)

#The best weight is "site"
residplot(m1.v1, newwd=F) 

#Homogeneity of variance improves , but still not great and there are some autocorrelation issues

#Applying log() to deal with heterogeneity and then check random structure

gls1 <- gls(f1b,
            method = "REML",
            weights = varIdent(form = ~ 1|type), #double check if this is still worth keeping
            data = soil1)

gls2 <- gls(f1b,
            method = "REML",
            data = soil1)

lme1 <- lme(f1b,
            random = ~ 1| site,
            method = "REML",
            data = soil1)

MuMIn::AICc(gls1, gls2, lme1) #log() considerably improves the models, lme1 has the best AICc

rm(f1a, f1b, gls1, gls2, lme2, m1.v1)

#Apply log() to the response variable####
soil1 <- soil1 %>% 
  mutate(Fe = log(Fe))

#Check residuals####
residplot(lme1, newwd = F) #homogeneity of variance issue with residuals, try varying weights
acf(residuals(lme1, type = "pearson")) #no autocorrelation

f1 <- Fe ~ treatment + type + position + treatment*type + type*position + treatment*position

lme2 <- lme(f1,
            random = ~ 1| site,
            weights = varIdent(form = ~1|site), 
            method = "REML",
            data = soil1)
lme3 <- lme(f1,
            random = ~ 1| site,
            weights = varIdent(form = ~1|type), 
            method = "REML",
            data = soil1)
lme4 <- lme(f1,
            random = ~ 1| site,
            weights = varIdent(form = ~1|treatment), 
            method = "REML",
            data = soil1)
lme5 <- lme(f1,
            random = ~ 1| site,
            weights = varIdent(form = ~1|position), 
            method = "REML",
            data = soil1)

AICc(lme1, lme2, lme3, lme4, lme5) # weight for site seems to improve model, check residuals

residplot(lme2, newwd = F) # much better
acf(residuals(lme2, type = "pearson")) # no autocorrelation


rm(lme1, lme2, lme3, lme4, lme5)

#Model selection ####
#Now, create model list for model selection, remember to run models by ML to compare fixed factors
f1 <- Fe ~ treatment + type + position + treatment*type + type*position + treatment*position
f2 <- Fe ~ treatment + type
f3 <- Fe ~ treatment + position
f4 <- Fe ~ position + type
f5 <- Fe ~ treatment + type + position
f6 <- Fe ~ treatment + type + position + treatment*type
f7 <- Fe ~ treatment + type + treatment*type
f8 <- Fe ~ treatment + type + position + type*position
f9 <- Fe ~ type + position + type*position
f10 <- Fe ~ treatment + type + position + treatment*position
f11 <- Fe ~ treatment + position + treatment*position
f12 <- Fe ~ treatment + type + position + treatment*type + type*position
f13 <- Fe ~ treatment + type + position + treatment*type + treatment*position
f14 <- Fe ~ treatment + type + position + treatment*position + treatment*type

m1 <- lme(f1,
          random = ~1|site,
          weights = varIdent(form = ~1|site),
          method = "ML",
          data = soil1)

m2 <- lme(f2,
          random = ~1|site,
          weights = varIdent(form = ~1|site),
          method = "ML",
          data = soil1)

m3 <- lme(f3,
          random = ~1|site,
          weights = varIdent(form = ~1|site),
          method = "ML",
          data = soil1)

m4 <-  lme(f4,
           random = ~1|site,
           weights = varIdent(form = ~1|site),
           method = "ML",
           data = soil1)

m5 <- lme(f5,
          random = ~1|site,
          weights = varIdent(form = ~1|site),
          method = "ML",
          data = soil1)

m6 <- lme(f6,
          random = ~1|site,
          weights = varIdent(form = ~1|site),
          method = "ML",
          data = soil1)

m7 <- lme(f7,
          random = ~1|site,
          weights = varIdent(form = ~1|site),
          method = "ML",
          data = soil1)

m8 <- lme(f8,
          random = ~1|site,
          weights = varIdent(form = ~1|site),
          method = "ML",
          data = soil1)

m9 <- lme(f9, 
          random = ~1|site,
          weights = varIdent(form = ~1|site),
          method = "ML",
          data = soil1)

m10 <- lme(f10,
           random = ~1|site,
           weights = varIdent(form = ~1|site),
           method = "ML",
           data = soil1)

m11 <- lme(f11,
           random = ~1|site,
           weights = varIdent(form = ~1|site),
           method = "ML",
           data = soil1)

m12 <- lme(f12,
           random = ~1|site,
           weights = varIdent(form = ~1|site),
           method = "ML",
           data = soil1)

m13 <- lme(f13,
           random = ~1|site,
           weights = varIdent(form = ~1|site),
           method = "ML",
           data = soil1)

m14 <- lme(f14,
           random = ~1|site,
           weights = varIdent(form = ~1|site),
           method = "ML",
           data = soil1)

#Run model selection
sel_iron <- model.sel(m1, m2, m3, m4, m5, m6, m7, m8, m9, m10, m11, m12, m13, m14)
#Look at the models in the table.
#Average the ones with delta <2, if just one, that is the best model
#m7 is the top model

rm(m1, m2, m3, m4, m5, m6, m7, m8, m9, m10, m11, m12, m13, m14)

#Final model ####
#Re-run by REML
m_iron <- lme(Fe ~ treatment + type + treatment*type,
              random = ~1|site,
              weights = varIdent(form = ~1|site),
              method = "REML",
              data = soil1)
summary(m_iron)

residplot(m_iron, newwd = F) #residuals OK
acf(residuals(m_iron)) # no autocorrelation


#Plot####
plot_data_fe <- m_iron %>% tidy_plus_plus(
                                          exponentiate = T, #exponentiate as the var was log transformed to deal with heterogeneity
                                          variable_labels = c("(Intercept)" = "Intercept",
                                                              # "position" = "Probe position",
                                                              "type" = "Carcass treatment",
                                                              "treatment" = "Exclusion treatment",
                                                              "type" = "Type"),
                                                      ) 
plot_data_fe <- plot_data_fe %>%
  filter(effect == "fixed") %>% 
  mutate(
    #A column to tell if the value was positive or negative. Remember here we have Post/pre ratio
    pos_neg = case_when(estimate < 1 & estimate >=0 ~ 1, # small positive value
                        estimate < 0 ~ -1, # negative value
                        estimate > 1 ~ 2,  # large positive value
                        estimate == 1 ~ 0) # reference value
  ) %>% 
  filter(term != "(Intercept)")

plot_iron <- plot_data_fe %>% 
  ggplot(mapping = aes(reorder(label, estimate), estimate,
                       col = pos_neg))+
  geom_point()+
  coord_flip()+
  geom_abline(intercept = 1,
              slope = 0,
              linetype = 2,
              colour = "black")+
  labs(x = "",
       y = "Estimate",
       title = "Model estimates for iron (Post/Baseline)")+
  geom_errorbar(aes(ymin = conf.low, ymax = conf.high), width = 0.2)+
  theme_classic()+
  theme(legend.position="none")

plot_iron

ggsave("figures/soil/models/m5_iron.png", plot_iron,
       width = 8,
       height = 6)

#Model results####
model_est <- m_iron %>% tidy_plus_plus(
  exponentiate = T, #exponentiate as the var was log transformed to deal with heterogeneity
  variable_labels = c("(Intercept)" = "Intercept",
                      #"position" = "Probe position",
                      "type" = "Carcass treatment",
                      "treatment" = "Exclusion treatment"
  )) %>% 
  dplyr::select(var_label, label, n_obs, estimate, std.error, conf.low, conf.high)%>% 
  mutate(across(.cols = c(3:6), ~round(.x, digits = 3))) %>% 
  rename("Variable" = "var_label",
         "Variable level" = "label",
         "Estimate" = "estimate", 
         "Number of observations" = "n_obs",
         "SE" = "std.error",
         "Lower CI (95%)" = "conf.low",
         "Upper CI (95%)" = "conf.high")

writexl::write_xlsx(model_est, "results/soil/m_iron.xlsx")

#Estimated model means####
model.rg <- update(ref_grid(m_iron), tran = "log") #to backtransform from log
emmeans_iron <- as.data.frame(emmeans(model.rg, ~ type | treatment, 
                                           type = "response"
))%>% 
  mutate(across(c(3,4,6,7), ~round(.x, digits = 1))) %>% 
  mutate(lower.SE = response - SE,
         upper.SE = response + SE)

writexl::write_xlsx(emmeans_iron, "results/soil/emmeans_m_iron.xlsx")

emmeans_iron_plot <- emmeans_iron %>% 
  ggplot(mapping = aes(factor(type,
                              levels = c("Control", "Single carcass", 'Mass mortality')),
                       response
  ))+
  geom_point()+
  labs(x = "",
       y = "Estimated ratio (Post/Baseline)",
       title = "Iron")+
  geom_errorbar(aes(ymin = lower.SE, ymax = upper.SE), width = 0.2)+
  theme_bw()+
  theme(legend.position="none")+
  facet_wrap(vars(treatment))+
  geom_text(label = emmeans_iron$response,
            nudge_x = 0.2,
            size = 3)

emmeans_iron_plot 

ggsave("figures/soil/models/emmeans_iron.png", emmeans_iron_plot,
       width = 9,
       height = 3)
