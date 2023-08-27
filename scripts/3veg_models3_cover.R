
#packages ####
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
library(emmeans)

# Load data ####

#please refer to previous script with data prep steps for general data exploration

veg <- read_excel("data/veg_analysis.xlsx")%>% 
  mutate(site = as.factor(site),
         type = as.factor(type),
         treatment = as.factor(treatment),
         transect = as.factor(transect),
         quadrat = as.factor(quadrat),
         treatment = gsub("Invertebrate exclusion", "Insect suppression", treatment)
           )

veg1 <- veg %>% 
  select(-forb, -grasses)

hist(veg1$cover)

#Model formulas ####
# f1 is the most complex model formula, therefore to be used for selecting optimal random structure
# Model formulas cannot be consistent with soil analysis formulas as the position variable is causing singularity issues in models
# I will use moisture as a variable in the models
# Moisture can affect plant growth, and differences in moisture may determine how much the plant grows

f1  <- cover ~ treatment + type + moisture + inorg_N + P + K + Mg + treatment*type + type*moisture + treatment*moisture
# Models for cover response ####
gls1 <- gls(f1,
            method = "REML",
            data = veg1
)


residplot(gls1, newwd=F) #homogeneity of variance a bit weird, bit of autocorrelation, resid normality a little off but not too bad


# First, try some variance structures to improve homogeneity of variance
gls1.v1 <- gls(f1,
               method = "REML",
               weights = varIdent(form = ~1| site), #change with different variables to see how they affect fit of model
               data = veg1
)
gls1.v2 <- gls(f1,
               method = "REML",
               weights = varIdent(form = ~1| type), #change with different variables to see how they affect fit of model
               data = veg1
)
gls1.v3 <- gls(f1,
               method = "REML",
               weights = varIdent(form = ~1| treatment), #change with different variables to see how they affect fit of model
               data = veg1
)

gls1.v4 <- gls(f1,
               method = "REML",
               weights = varIdent(form = ~1| transect), #change with different variables to see how they affect fit of model
               data = veg1
)
gls1.v5 <- gls(f1,
               method = "REML",
               weights = varIdent(form = ~1| quadrat), #change with different variables to see how they affect fit of model
               data = veg1
)

MuMIn::AICc(gls1, gls1.v1, gls1.v2, gls1.v3, gls1.v4, gls1.v5) #v2 = Type is the best, remove all others
rm(gls1, gls1.v1, gls1.v3, gls1.v4, gls1.v5)

#validation plot
residplot(gls1.v2, newwd=F)  #little improvement

# Adding random effects###
# Now, see if adding a random effect improves models

lme1 <- lme(f1,
            random = ~1| site/transect/quadrat ,
            method = "REML",
            weights = varIdent(form = ~1| type),
            data = veg1)


# Let's try with another different random structure
lme2 <- lme(f1,
            random = ~1| transect/quadrat ,
            method = "REML",
            weights = varIdent(form = ~1| type),
            data = veg1)

lme3 <- lme(f1,
            random = ~1| quadrat ,
            method = "REML",
            weights = varIdent(form = ~1| type),
            data = veg1)

#Check if leaving the weights is a good idea
lme4 <- lme(f1,
            random = ~1| site/transect/quadrat ,
            method = "REML",
            data = veg1)

MuMIn::AICc(gls1.v2, lme1, lme2, lme3, lme4) #lme1 has the best fit

#Check residuals and autocorrelation
residplot(lme1, newwd=F, level = 3)  
acf(residuals(lme1, type = "pearson")) #improved a lot

#see if adding a correlation structure helps
lme5 <- lme(f1,
            random = ~1| site/transect/quadrat ,
            method = "REML",
            weights = varIdent(form = ~1| site),
            correlation = corCompSymm(form = ~1|site/transect/quadrat),
            data = veg1)

AICc(lme1, lme5) #No

#Keep lme1 as optimal structure
rm(lme1, lme2, lme3, lme4, lme5)


#Model selection, fixed factors (by ML) ####
f1  <- cover ~ treatment + type + moisture + inorg_N + P + K + Mg + treatment*type + type*moisture + treatment*moisture

m1  <- lme(f1,
           random = ~1|site/transect/quadrat,
           weights = varIdent(form = ~1|site),
           correlation = corCompSymm(form = ~1|site/transect/quadrat),
           method = "ML",
           data = veg1)

summary(m1)

#Drop one variable, and then choose a model with best AIC
m1.dropped <- drop1(m1)
m1.dropped

#Drop Potassium
m2 <- lme(cover ~ treatment + type + moisture + inorg_N + P + Mg + treatment*type + type*moisture + treatment*moisture,
          random = ~1|site/transect/quadrat,
          weights = varIdent(form = ~1|site),
          correlation = corCompSymm(form = ~1|site/transect/quadrat),
          method = "ML",
          data = veg1)
summary(m2)

m2.dropped <- drop1(m2)
m2.dropped #doesn't look like dropping any other variables improves the model

#m2 might be the best fit

rm(m1, m2, m1.dropped, m2.dropped)

#Final model####
#Re-run by "REML"

m_cover  <- lme(cover ~ treatment + type + moisture + inorg_N + P + Mg + treatment*type + type*moisture + treatment*moisture,
               random = ~1|site/transect/quadrat,
               weights = varIdent(form = ~1|site),
               correlation = corCompSymm(form = ~1|site/transect/quadrat),
               method = "REML",
               data = veg1)
summary(m_cover)

residplot(m_cover, newwd = F) # residuals OK
acf(residuals(m_cover, type = "pearson")) #good

#Plot####
plot_data_cover <- m_cover %>% tidy_plus_plus(
  tidy_fun = broom.helpers::tidy_parameters,
  effects = "fixed",
  add_reference_rows = T,
  add_estimate_to_reference_rows = T,
  exponentiate = F, #data was not log transformed
  variable_labels = c(#"(Intercept)" = "Intercept",
    "moisture" = "Moisture",
    "type" = "Carcass treatment",
    "treatment" = "Exclusion treatment",
    "type" = "Type",
    "inorg_N" = "Inorganic nitrogen",
    "P" = "Phosphorus",
    "Mg" = "Magnesium")
  ) %>% 
  select(var_label, label,  estimate, std.error, conf.low, conf.high, statistic, p.value)%>% 
  mutate(across(3:8, ~round(.x, digits = 2))) %>% 
  mutate(label = str_replace_all(label, "treatmentOpen", "Open"),
         label = str_replace_all(label, "treatmentVertebrate exclusion", "Vertebrate exclusion")) %>% 
  add_row(var_label = "Exclusion treatment", label = "Insect suppression", estimate = 0, 
          std.error = NA, conf.low = NA, conf.high = NA, statistic = NA, p.value = NA,
          .after = 2)

plot_data_cover <- plot_data_cover %>%
  mutate(
    #A column to tell if the value was positive or negative. Remember here we have Post/pre ratio
    pos_neg = case_when(estimate > 0 ~ 1, # positive value
                        estimate < 0 ~ -1, # negative value
                        estimate == 0 ~ 0) # reference value
  ) 

plot_cover <- plot_data_cover %>% 
  ggplot(mapping = aes(reorder(label, estimate), estimate,
                       col = pos_neg))+
  geom_point()+
  coord_flip()+
  geom_abline(intercept = 0, #reference line
              slope = 0,
              linetype = 2,
              colour = "black")+
  labs(x = "",
       y = "Estimate",
       title = "")+
  geom_errorbar(aes(ymin = conf.low, ymax = conf.high), width = 0.2)+
  theme_classic()+
  theme(legend.position="none")

plot_cover

ggsave("figures/veg/models/m3_cover.png", plot_cover, 
       height = 6,
       width = 8)


#Model results####
model_est <- m_cover %>% tidy_plus_plus(
  tidy_fun = broom.helpers::tidy_parameters,
  effects = "fixed",
  add_reference_rows = T,
  add_estimate_to_reference_rows = T,
  exponentiate = F, 
  variable_labels = c("moisture" = "Moisture",
                      "type" = "Carcass treatment",
                      "treatment" = "Exclusion treatment",
                      "type" = "Type",
                      "inorg_N" = "Inorganic nitrogen",
                      "P" = "Phosphorus",
                      "Mg" = "Magnesium")) %>% 
  dplyr::select(var_label, label, n_obs, estimate, std.error, conf.low, conf.high)%>% 
  mutate(t.value = estimate / std.error,
         across(.cols = c(3:8), ~round(.x, digits = 3)),
         .after = std.error) %>% 
  rename("Variable" = "var_label",
         "Variable level" = "label",
         "Estimate" = "estimate", 
         "Number of observations" = "n_obs",
         "SE" = "std.error",
         "Lower CI (95%)" = "conf.low",
         "Upper CI (95%)" = "conf.high")

writexl::write_xlsx(model_est, "results/veg/m_cover.xlsx")

#Estimated model means####
emmeans_cover <- as.data.frame(emmeans(m_cover, ~ type | treatment | inorg_N | P | Mg, 
                                       type = "response")) %>% 
  mutate(across(c(3:5), ~round(.x, digits = 2))) %>% 
  mutate(across(c(6, 7, 9, 10), ~round(.x*100, digits = 1))) %>% 
  mutate(lower.SE = emmean - SE,
         upper.SE = emmean + SE)
writexl::write_xlsx(emmeans_cover, "results/veg/emmeans_m_cover.xlsx")


emmeans_cover_plot <- emmeans_cover %>% 
  ggplot(mapping = aes(factor(type,
                              levels = c("Control", "Single carcass", 'Mass mortality')),
                       emmean)
         )+
  geom_point()+
  labs(x = "",
       y = "Estimated change (Post - Baseline)",
       title = "Overall vegetation cover")+
  geom_errorbar(aes(ymin = emmean+SE, ymax = emmean-SE), width = 0.2)+
  geom_hline(yintercept = 0, linetype = "dotted", color = "red", size = 0.5) +
  theme_classic()+
  theme(legend.position="none",
        panel.border = element_rect(fill = NA, colour = "black", linewidth = 1))+
  facet_wrap(vars(treatment)) +
  # scale_x_discrete(labels = c("Open","Invertebrate ex",  "Vertebrate ex"))+
  geom_text(label = paste(emmeans_cover$emmean, "%", sep = ""),
            nudge_x = 0.3,
            size = 3.5)

emmeans_cover_plot


ggsave("figures/veg/models/emmeans_cover.png", emmeans_cover_plot,
       width = 9,
       height = 3)


veg_plots <- ggarrange(emmeans_cover_plot, emmeans_grasses_plot, emmeans_forb_plot,
                       ncol = 1,
                       nrow = 3,
                       labels = "AUTO")

ggsave("figures/emmeans_veg.png", veg_plots,
       width = 9,
       height = 9)