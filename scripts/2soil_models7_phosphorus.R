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
  select(plot_id, sample_n, site, type, treatment, position, P) %>% 
  mutate(site = as.factor(site),
         type = as.factor(type),
         treatment = as.factor(treatment),
         position = as.factor(position)
  )
#See the data prep script within this project for all the data exploration

#phosphorus ####

# Formulas. f1a, f1b, f1 are the formulas including all terms possible. 
# These are used to find the optimal variance structure, then to find the optimal random structure
# The rest of the formulas will be updated later when the above has been accomplished

f1a <- P ~ treatment + type + position + treatment*type + type*position + treatment*position

#Log transform data
f1b <- log(P) ~ treatment + type + position + treatment*type + type*position + treatment*position

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
rm(m1, m1.v2, m1.v1, m1.v4)

#The best weight is "treatment"
residplot(m1.v3, group = "type", newwd=F) 

#Homogeneity of variance improves , but still not great and there are some autocorrelation issues

#Applying log() to deal with heterogeneity and then check random structure

gls1 <- gls(f1b,
            method = "REML",
            weights = varIdent(form = ~ 1|treatment), #double check if this is still worth keeping
            data = soil1)

gls2 <- gls(f1b,
            method = "REML",
            data = soil1)

lme1 <- lme(f1b,
            random = ~ 1| site,
            method = "REML",
            data = soil1)

lme2 <- lme(f1b,
            random = ~ 1| site,
            weights = varIdent(form = ~1|treatment), #double check if this is still worth keeping
            method = "REML",
            data = soil1)

MuMIn::AICc(gls1, gls2, lme1, lme2) #log() considerably improves the models, gls2 and lme1 have the best AICc

rm(f1a, f1b, gls1, gls2, lme2, m1.v1)

#Apply log() to the response variable####
soil1 <- soil1 %>% 
  mutate(P = log(P))

#check residuals####
residplot(lme1, newwd = F) #residuals look pretty good
acf(residuals(lme1, type = "pearson")) #no autocorrelation

f1 <- P ~ treatment + type + position + treatment*type + type*position + treatment*position


gls2 <- gls(f1,
            method = "REML",
            data = soil1)  
residplot(gls2, newwd = F) #residuals look pretty good, but there's some autocorrelation

#therefore lme1 is better
#lme1 also has good homogeneity of variance, so no need for weights

rm(lme1, gls2)

#Model selection ####
#Now, create model list for model selection, remember to run models by ML to compare fixed factors
f1 <- P ~ treatment + type + position + treatment*type + type*position + treatment*position
f2 <- P ~ treatment + type
f3 <- P ~ treatment + position
f4 <- P ~ position + type
f5 <- P ~ treatment + type + position
f6 <- P ~ treatment + type + position + treatment*type
f7 <- P ~ treatment + type + treatment*type
f8 <- P ~ treatment + type + position + type*position
f9 <- P ~ type + position + type*position
f10 <- P ~ treatment + type + position + treatment*position
f11 <- P ~ treatment + position + treatment*position
f12 <- P ~ treatment + type + position + treatment*type + type*position
f13 <- P ~ treatment + type + position + treatment*type + treatment*position
f14 <- P ~ treatment + type + position + treatment*position + treatment*type

m1 <- lme(f1,
          random = ~1|site,
          method = "ML",
          data = soil1)

m2 <- lme(f2,
          random = ~1|site,
          method = "ML",
          data = soil1)

m3 <- lme(f3,
          random = ~1|site,
          method = "ML",
          data = soil1)

m4 <-  lme(f4,
           random = ~1|site,
           method = "ML",
           data = soil1)

m5 <- lme(f5,
          random = ~1|site,
          method = "ML",
          data = soil1)

m6 <- lme(f6,
          random = ~1|site,
          method = "ML",
          data = soil1)

m7 <- lme(f7,
          random = ~1|site,
          method = "ML",
          data = soil1)

m8 <- lme(f8,
          random = ~1|site,
          method = "ML",
          data = soil1)

m9 <- lme(f9, 
          random = ~1|site,
          method = "ML",
          data = soil1)

m10 <- lme(f10,
           random = ~1|site,
           method = "ML",
           data = soil1)

m11 <- lme(f11,
           random = ~1|site,
           method = "ML",
           data = soil1)

m12 <- lme(f12,
           random = ~1|site,
           method = "ML",
           data = soil1)

m13 <- lme(f13,
           random = ~1|site,
           method = "ML",
           data = soil1)

m14 <- lme(f14,
           random = ~1|site,
           method = "ML",
           data = soil1)

#Run model selection
sel_phosphorus <- model.sel(m1, m2, m3, m4, m5, m6, m7, m8, m9, m10, m11, m12, m13, m14)
#lookat the models in the table.
#Average the ones with delta <2, if just one, that is the best model
#m4, m3 are the best models

m_avgd <- model.avg(sel_phosphorus, subset = delta <2)
summary(m_avgd)

rm(m1, m2, m3, m4, m5, m6, m7, m8, m9, m10, m11, m12, m13, m14, m_avgd)

#Final model ####
#Re-run by REML
m_phosphorus <- lme(P ~ position + treatment + type,
                   random = ~1|site,
                   method = "REML",
                   data = soil1)
summary(m_phosphorus)

residplot(m_phosphorus, newwd = F) #residuals OK
acf(residuals(m_phosphorus)) # no autocorrelation


#Plot####

plot_data_p <- m_phosphorus %>% tidy_plus_plus(
                                      exponentiate = T, #exponentiate as the var was log transformed to deal with heterogeneity
                                      variable_labels = c("(Intercept)" = "Intercept",
                                      "position" = "Probe position",
                                      "type" = "Carcass treatment",
                                      "treatment" = "Exclusion treatment",
                                      "type" = "Type")
                                    ) 
plot_data_p <- plot_data_p %>%
  filter(effect == "fixed") %>% 
  mutate(
    #A column to tell if the value was positive or negative. Remember here we have Post/pre ratio
    pos_neg = case_when(estimate < 1 & estimate >=0 ~ 1, # small positive value
                        estimate < 0 ~ -1, # negative value
                        estimate > 1 ~ 2,  # large positive value
                        estimate == 1 ~ 0) # reference value
  ) %>% 
  filter(term != "(Intercept)")

plot_phosphorus <- plot_data_p %>% 
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
       title = "Model estimates for phosphorus (Post/Baseline)")+
  geom_errorbar(aes(ymin = conf.low, ymax = conf.high), width = 0.2)+
  theme_classic()+
  theme(legend.position="none")

plot_phosphorus

ggsave("figures/soil/models/m7_phosphorus.png", plot_phosphorus,
       width = 8,
       height = 6)
#Model results####
model_est <- m_phosphorus %>% tidy_plus_plus(
  exponentiate = T, #exponentiate as the var was log transformed to deal with heterogeneity
  variable_labels = c("(Intercept)" = "Intercept",
                      "position" = "Probe position",
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

writexl::write_xlsx(model_est, "results/soil/m_phosphorus.xlsx")


#Estimated model means####
model.rg <- update(ref_grid(m_phosphorus), tran = "log") #to backtransform from log
emmeans_phosphorus <- as.data.frame(emmeans(model.rg, ~ position | type | treatment, 
                                            type = "response"
                                            ))%>% 
  mutate(across(c(4,5,7,8), ~round(.x, digits = 1))) %>% 
  mutate(lower.SE = response - SE,
         upper.SE = response + SE)

writexl::write_xlsx(emmeans_phosphorus, "results/soil/emmeans_m_phosphorus.xlsx")

emmeans_phosphorus_plot <- emmeans_phosphorus %>% 
  ggplot(mapping = aes(factor(type,
                              levels = c("Control", "Single carcass", 'Mass mortality')), response
  ))+
  geom_point()+
  labs(x = "",
       y = "Estimated ratio (Post/Baseline)",
       title = "Phosphorus")+
  geom_errorbar(aes(ymin = response -SE, ymax = response + SE), width = 0.2)+
  theme_bw()+
  theme(legend.position="none",
        axis.text.x = element_text(size = 7))+
  facet_grid(position~treatment)+
  geom_text(label = emmeans_phosphorus$response,
            nudge_x = 0.2,
            size = 3)

emmeans_phosphorus_plot 

ggsave("figures/soil/models/emmeans_phosphorus.png", emmeans_phosphorus_plot,
       width = 8,
       height = 4)

plot_pmgmn <- ggarrange(emmeans_phosphorus_plot, 
                        emmeans_magnesium_plot,
                        emmeans_manganese_plot,
                        labels = "AUTO",
                        ncol = 1,
                        nrow = 3)

ggsave("figures/emmeans_pmgmn.png", plot_pmgmn,
       width = 8,
       height = 8)