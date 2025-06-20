# RQ1: In what ways is the aquaculture land transition related to salinity prevalence and history of storm surges,
# accounting for place-invariant factors such as distance/elevation from the sea or soil conditions, or time based shocks such as policy changes, inflation, or market conditions

# RQ2: How do these relationships vary by states? 

# Optional RQ3: Is there a predictive lag between surge impact and land change outcome? 


# Aquaculture Analysis: Panel Data Regression Models
# This script analyzes the relationship between aquaculture land change, storm surges, salinity, and other socio-environmental factors


# ---------------------------------------
# DATA AND PREP
# ---------------------------------------
# Load required packages
library(plm)       # For panel data models
library(ggplot2)   # For visualization
library(dplyr)     # For data manipulation
library(stargazer) # For regression table output
library(lmtest)    # For statistical tests
library(sandwich)  # For robust standard errors
library(tidyr)     # Fill missing
library(tidyverse)
library(modelsummary) # for summary tables
library(fixest)    # for fixed effects with feols()
library(zoo) #for rollmean() for calculating smoother percentage changes over 3-5 years
library(slider)getwd()
library(margins) # For visualizing 



data <- read.csv(unz("outputs/aqua_salinity_surge_1990-2025.zip", 
                           "aqua_salinity_surge_1990-2025.csv"))

head(data)
summary(data)

# Check for missing values
missing_values <- colSums(is.na(data))
print(missing_values)
# Many pct_change_aqua and pct_change_salinity are NA. Might not want to use that for main analysis.  
# Convert pct change to categorical and then use in analysis

# Make year numeric for all variable creation and ordering
data$Year <- as.numeric(as.character(data$Year))

# Convert year to factor for fixed effects
data$year_factor <- as.factor(data$Year)

summary(data$avg_salinity_5yr)

# Set TN as a reference state
data$State <- relevel(factor(data$State), ref = "TN")


# ---------------------------------------
# REGRESSION MODELS : BASIC RELATIONSHIPS
# ---------------------------------------
# Two-way fixed effect models 

# SET A : Without and with controls 
# 1. Y = Salinity, X = surge (without controls except time and unit, with controls)
regA1_1 <- feols(avg_salinity_5yr ~ postSurge | Year, 
                data = data, vcov = ~UniqueID)
summary(regA1_1)

regA1_2 <- feols(avg_salinity_5yr ~ postSurge | UniqueID + Year, 
                data = data, vcov = ~UniqueID)
summary(regA1_2)

# When comparing across places, Villages in post-surge years have, on average, 0.608 units higher 5-year average salinity, compared to others, controlling only for time trends. 
# When comparing villages to themselves overtime, the relationship of salinity with surges is still positive and significant, although smaller. 
# This suggests that while places that are affected by surges have a higher baseline salinity (i.e. are more likely to be saline), but exprience relatively more salinity after being affected by storms. 


# 2. Y = Aquaculture, X = Surge (without controls except time and unit, with controls)
regA2_1 <- feols(Aqua_perc ~ postSurge | Year, 
                 data = data, vcov = ~UniqueID)
summary(regA2_1)

regA2_2 <- feols(Aqua_perc ~ postSurge | UniqueID + Year, 
                 data = data, vcov = ~UniqueID)
summary(regA2_2)

# When comparing across places, places affected by surges are highly likely to have higher aquaculture. 
# but when comparing within places, Within the same village, aquaculture increases by ~0.16 percentage points after a surge.
# This is accounting for place-based time-invariant village characteristics and time shocks that affect all villages (policies, inflation, prices, etc.), and therefore likely a causal effect.  


# 3. Y = Aquaculture, X = Salinity in previous period (without controls except time and unit, with controls)
regA3_1 <- feols(Aqua_perc ~ avg_salinity_5yr | Year, 
                 data = data, vcov = ~UniqueID)
summary(regA3_1)

regA3_2 <- feols(Aqua_perc ~ avg_salinity_5yr | UniqueID + Year, 
                 data = data, vcov = ~UniqueID)
summary(regA3_2)



# 4. Y = Aquaculture, X = Aqua time - 1 (without controls except time and unit, with controls)
regA4_1 <- feols(Aqua_perc ~ Lag_Aqua | Year, 
                 data = data, vcov = ~UniqueID)
summary(regA4_1)

regA4_2 <- feols(Aqua_perc ~ Lag_Aqua | UniqueID + Year, 
                 data = data, vcov = ~UniqueID)
summary(regA4_2)

# Aquaculture in the previous period remains consistently high and significant predictor or aquaculture in the current period



# 5. Y = Aquaculture, X = Surge, X = Salinity, X = Aqua-1 (without controls except time and unit, with controls)
regA5_1 <- feols(Aqua_perc ~ postSurge + avg_salinity_5yr + Lag_Aqua | Year, 
                 data = data, vcov = ~UniqueID)
summary(regA5_1)

regA5_2 <- feols(Aqua_perc ~ postSurge + avg_salinity_5yr + Lag_Aqua | UniqueID + Year, 
                 data = data, vcov = ~UniqueID)
summary(regA5_2)

# When including all predictors, they all seem to have a strong and significant relationship with aquaculture, both within and across comparisons. 


models <- list()
models[['Salinity (Time FE)']] <- regA1_1
models[['Salinity (Village+Time FE)']] <- regA1_2

models[['AC (Time FE)']] <- regA2_1
models[['AC (Village+Time FE)']] <- regA2_2

models[['AC (time t) (only Time FE)']] <- regA3_1
models[['AC (time t) (Village+Time FE)']] <- regA3_2

models[['AC (in time t) (Time FE)']] <- regA4_1
models[['AC (in time t) (Village+Time FE)']] <- regA4_2

models[['AC (only Time FE)']] <- regA5_1
models[['AC (All controls+Village+Time FE)']] <- regA5_2

msummary(models, stars = c('*' = .1, '**' = .05, '***' = .01),gof_omit=c("BIC|AIC|RMSE|R2 Within Adj."),coef_omit=c("(Intercept)"), filename = 'table.rtf')




# Set B : With interactions 
# 1. Y = Aquaculture, X = Surge, X = Salinity, X = Aqua-1 

regB1_1 <- feols(Aqua_perc ~ postSurge + avg_salinity_5yr + Lag_Aqua + postSurge * avg_salinity_5yr | Year, 
                 data = data, vcov = ~UniqueID)
summary(regB1_1)

regB1_2 <- feols(Aqua_perc ~ postSurge + avg_salinity_5yr + Lag_Aqua + postSurge * avg_salinity_5yr | UniqueID + Year, 
                 data = data, vcov = ~UniqueID)
summary(regB1_2)


models <- list()
models[['AC (All controls+Time FE)']] <- regB1_1
models[['AC (All controls+Village+Time FE)']] <- regB1_2

msummary(models, stars = c('*' = .1, '**' = .05, '***' = .01),gof_omit=c("BIC|AIC|RMSE|R2 Within Adj."),coef_omit=c("(Intercept)"), filename = 'table.rtf')




# Set C : State variation

# 1. Y = Salinity, X = surge (without controls except time and unit, with controls)
regC1_1 <- feols(avg_salinity_5yr ~ postSurge + postSurge * State | Year, 
                 data = data, vcov = ~UniqueID)
summary(regC1_1)

regC1_2 <- feols(avg_salinity_5yr ~ postSurge + postSurge * State | UniqueID + Year, 
                 data = data, vcov = ~UniqueID)
summary(regC1_2)


# 2. Y = Aquaculture, X = Surge (without controls except time and unit, with controls)
regC2_1 <- feols(Aqua_perc ~ postSurge + postSurge * State | Year, 
                 data = data, vcov = ~UniqueID)
summary(regC2_1)

regC2_2 <- feols(Aqua_perc ~ postSurge + postSurge * State | UniqueID + Year, 
                 data = data, vcov = ~UniqueID)
summary(regC2_2)



# 3. Y = Aquaculture, X = Salinity in previous period (without controls except time and unit, with controls)
regC3_1 <- feols(Aqua_perc ~ avg_salinity_5yr + avg_salinity_5yr * State | Year, 
                 data = data, vcov = ~UniqueID)
summary(regC3_1)

regC3_2 <- feols(Aqua_perc ~ avg_salinity_5yr + avg_salinity_5yr * State | UniqueID + Year, 
                 data = data, vcov = ~UniqueID)
summary(regC3_2)



# 4. Y = Aquaculture, X = Aqua time - 1 (without controls except time and unit, with controls)
regC4_1 <- feols(Aqua_perc ~ Lag_Aqua + Lag_Aqua * State | Year, 
                 data = data, vcov = ~UniqueID)
summary(regC4_1)

regC4_2 <- feols(Aqua_perc ~ Lag_Aqua + Lag_Aqua * State | UniqueID + Year, 
                 data = data, vcov = ~UniqueID)
summary(regC4_2)



models <- list()
models[['Salinity (Time FE)']] <- regC1_1
models[['Salinity (Village+Time FE)']] <- regC1_2

models[['AC (Time FE)']] <- regC2_1
models[['AC (Village+Time FE)']] <- regC2_2

models[['AC (time t) (only Time FE)']] <- regC3_1
models[['AC (time t) (Village+Time FE)']] <- regC3_2

models[['AC (in time t) (Time FE)']] <- regC4_1
models[['AC (in time t) (Village+Time FE)']] <- regC4_2

msummary(
  models,
  stars = c('*' = .1, '**' = .05, '***' = .01),
  gof_omit = c("BIC|AIC|RMSE|R2 Within Adj."),
  coef_omit = c("(Intercept)"),
  coef_map = c(
    "postSurge" = "Post-Surge (TN)",
    "postSurge:StateAP" = "Post-surge (AP wrt TN)", 
    "postSurge:StateOD" = "Post-surge (OD wrt TN)",
    "avg_salinity_5yr" = "Persisent Salinity (TN)",
    "avg_salinity_5yr:StateAP" = "Persisent Salinity (AP wrt TN)",
    "avg_salinity_5yr:StateOD" = "Persisent Salinity (OD wrt TN)",
    "Lag_Aqua" = "AC in t-1 (TN)",
    "Lag_Aqua:StateAP" = "AC in t-1 (AP wrt TN)",
    "Lag_Aqua:StateOD" = "AC in t-1 (OD wrt TN)",
    "StateAP" = "Baseline AP (wrt TN)",
    "StateOD" = "Baseline OD (wrt TN)"
  ),
  filename = 'table.rtf'
)


# Set D : State variations and Interactions

# 1. Y = Aquaculture, X = Surge, X = Salinity, X = Aqua-1 (with one interaction at a time)
regD1_1 <- feols(Aqua_perc ~ postSurge + avg_salinity_5yr + Lag_Aqua + postSurge * State | Year, 
                 data = data, vcov = ~UniqueID)
summary(regD1_1)

regD1_2 <- feols(Aqua_perc ~ postSurge + avg_salinity_5yr + Lag_Aqua + postSurge * State | UniqueID + Year, 
                 data = data, vcov = ~UniqueID)
summary(regD1_2)

# There is a strong heterogeneity in surge responses across states. Storm surges strongly increase aquaculture in AP, weakly decrease it in OD, and modestly increase it in TN.
# Once village-specific baselines are controlled, the surge effect in TN is strongest, smaller in AP, but extremely low in OD (0.012).
# Cross-state differences persist, but magnitude and direction flip compared to model without FE.



regD1_3 <- feols(Aqua_perc ~ postSurge + avg_salinity_5yr + Lag_Aqua + avg_salinity_5yr * State | Year, 
                 data = data, vcov = ~UniqueID)
summary(regD1_3)

regD1_4 <- feols(Aqua_perc ~ postSurge + avg_salinity_5yr + Lag_Aqua + avg_salinity_5yr * State | UniqueID + Year, 
                 data = data, vcov = ~UniqueID)
summary(regD1_4)

# Comparing places across the data, in all states higher salinity is related with higher aquaculture, although AP seems to be most responsive to salinity and OD the least. 
# Once unobserved village-level time invarying factors (such as dist or elevation from the sea), etc are controlled for, TN shows a small negative effect of salinity, whereas AP and OD show positive association. 


models <- list()
models[['AC (surge by states) (Time FE)']] <- regD1_1
models[['AC (surge by states) (Village+Time FE)']] <- regD1_2
models[['AC (salinity by states) (Time FE)']] <- regD1_3
models[['AC (salinity by states) (Village+Time FE)']] <- regD1_4

msummary(
  models,
  stars = c('*' = .1, '**' = .05, '***' = .01),
  gof_omit = c("BIC|AIC|RMSE|R2 Within Adj."),
  coef_omit = c("(Intercept)"),
  coef_map = c(
    "postSurge" = "Post-Surge (TN or all)",
    "postSurge:StateAP" = "Post-surge (AP wrt TN)", 
    "postSurge:StateOD" = "Post-surge (OD wrt TN)",
    "avg_salinity_5yr" = "Persisent Salinity (TN or all)",
    "avg_salinity_5yr:StateAP" = "Persisent Salinity (AP wrt TN)",
    "avg_salinity_5yr:StateOD" = "Persisent Salinity (OD wrt TN)",
    "Lag_Aqua" = "AC in t-1 (TN)",
    "Lag_Aqua:StateAP" = "AC in t-1 (AP wrt TN)",
    "Lag_Aqua:StateOD" = "AC in t-1 (OD wrt TN)",
    "StateAP" = "Baseline AP (wrt TN)",
    "StateOD" = "Baseline OD (wrt TN)"
  ),
  filename = 'table.rtf'
)

# TN aquaculture development seems to be most responsive to shocks, but less so to salinity stress. 
# AP seems to be strategically adaptive to long-term salinity. 
# OD shows least responsiveness to either shocks or stresses, likely due to infrastructure, policy or economic constraints. 




# Set E : Three-way interactions and non-linearity 
# Try triple interaction: 
regE1_1 <- feols(Aqua_perc ~ postSurge + avg_salinity_5yr + Lag_Aqua + postSurge * avg_salinity_5yr * State | UniqueID + Year, 
                 data = data, vcov = ~UniqueID)
summary(regE1_1)

# Interpretation: 
# TN - Surges have a big positive effect on aquaculture; long-run salinity has a small negative effect, salinity may amplify surge effect, but not significantly. 
# AP - Surges alone are less impacxtful as compared to TN, but salinity is strongly positively associated with aquaculture. Surges in high saline conditions tend to trigger aquaculture expansion the most in AP. 
# OD - Surges have a weak or even negative effect on aqua expansion. Salinity has a more positive response to aquaculture compared to TN. Surges in high saline places have small negative effect, however not significant. 


# Test for non-linear relationship with salinity 
data_clean <- data %>% filter(!is.na(avg_salinity_5yr))
regE1_2 <- feols(Aqua_perc ~ postSurge + poly(avg_salinity_5yr, 2) + Lag_Aqua + postSurge * poly(avg_salinity_5yr, 2) * State | UniqueID + Year, 
                 data = data_clean, vcov = ~UniqueID)
summary(regE1_2)

# Linear effect not significant alone
# Strong curvature effect — suggests increasing returns to salinity for aquaculture in base state.
# Some state-specific nonlinear salinity effects — e.g., salinity effects are stronger in AP and OD in higher-order terms.
# A massive, positive coefficient — suggests that in AP, aquaculture after a surge responds very strongly and nonlinearly to salinity - maybe visualise it? 

# 5-year average salinity has a nonlinear (convex) effect, meaning aquaculture increases at higher levels of salinity.
# In AP - Aquaculture is highly sensitive to salinity after a surge (strong positive interaction), implying surge and salinity may be joint enablers of expansion.
# In OD - Surge × salinity interactions are not significant. Salinity increases aquaculture moderately, but surge may disrupt or delay response.
# Lagged aquaculture is very predictive — suggesting strong path dependence or persistence in village-level aquaculture.

models <- list()
models[['AC (three-way interaction) (All controls+Two way FE)']] <- regE1_1
models[['AC (non-linearity) (All controls+Two way FE)']] <- regE1_2

msummary(models, stars = c('*' = .1, '**' = .05, '***' = .01),gof_omit=c("BIC|AIC|RMSE|R2 Within Adj."),coef_omit=c("(Intercept)"), filename = 'table.rtf')




# ------------------------------
# Y = Salinity, X = Aqua in t-1
# ------------------------------
# Calculate avg salinity for t-1 to t-5 (Note: current variable avg_salinity_5yr is t to t-4)
data <- data %>%
  mutate(Year = as.numeric(Year)) %>%
  arrange(UniqueID, Year) %>%
  group_by(UniqueID) %>%
  complete(Year = full_seq(Year, 1)) %>%
  mutate(
    Saline_perc_norm_lagged = lag(Saline_perc_norm, 1),  # shift by one year
    Saline_past5yr = rollapply(Saline_perc_norm_lagged, 
                                                width = 5, 
                                                FUN = mean, 
                                                fill = NA, 
                                                align = "right", 
                                                na.rm = TRUE, 
                                                partial = TRUE)
  ) %>%
  ungroup()



# Calculate future 5 year avg salinity 
data <- data %>%
  arrange(UniqueID, Year) %>%
  group_by(UniqueID) %>%
  mutate(
    # Shift salinity forward by 1 year so current Aqua matches future Saline
    Saline_lead = dplyr::lead(Saline_perc_norm, 1),
    
    # Calculate future 5-year average salinity: t+1 to t+5
    future_avg_salinity = rollapply(
      Saline_perc_norm,
      width = 5,
      align = "left",  # this ensures window starts at t
      FUN = function(x) mean(x[-1], na.rm = TRUE),  # exclude current year
      fill = NA,
      partial = TRUE
    )
  ) %>%
  ungroup()


# Restrict to villages only with aqua 
villages_with_aqua <- data %>%
  group_by(UniqueID) %>%
  summarize(has_aqua = any(!is.na(Aqua_perc) & Aqua_perc > 0)) %>%
  filter(has_aqua) %>%
  pull(UniqueID)

data_filtered <- data %>%
  filter(UniqueID %in% villages_with_aqua)


# 1. Y = Current Salinity, X = Lag Aqua (without controls except time and unit, with controls)
regF1_1 <- feols(Saline_perc_norm  ~ Lag_Aqua | UniqueID + Year, 
                 data = data_filtered, vcov = ~UniqueID)
summary(regF1_1)

regF1_2 <- feols(Saline_perc_norm  ~ Lag_Aqua + postSurge | UniqueID + Year, 
                 data = data_filtered, vcov = ~UniqueID)
summary(regF1_2)

# Adding postSurge has no impact on the model or the coefficient of past aquaculture 

models <- list()
models[['Salinity (Village+Time FE)']] <- regF1_1
models[['Salinity (All controls+Village+Time FE)']] <- regF1_2
msummary(models, stars = c('*' = .1, '**' = .05, '***' = .01),gof_omit=c("BIC|AIC|RMSE|R2 Within Adj."),coef_omit=c("(Intercept)"), filename = 'table.rtf')




# 2. Y = Future avg salinity, X = Current Aqua
regF2_1 <- feols(future_avg_salinity ~ Aqua_perc | UniqueID + Year, 
                 vcov = ~UniqueID, data = data_filtered)

summary(regF2_1)

regF2_2 <- feols(future_avg_salinity ~ Aqua_perc + postSurge | UniqueID + Year, 
                 vcov = ~UniqueID, data = data_filtered)

summary(regF2_2)

models <- list()
models[['Future Salinity (Village+Time FE)']] <- regF2_1
models[['Future Salinity (All controls+Village+Time FE)']] <- regF2_2
msummary(models, stars = c('*' = .1, '**' = .05, '***' = .01),gof_omit=c("BIC|AIC|RMSE|R2 Within Adj."),coef_omit=c("(Intercept)"), filename = 'table.rtf')



# These results provide strong empirical evidence of a lagged ecological feedback from human adaptation (aquaculture expansion) to system degradation (salinization). Even after removing the "zero-aquaculture" villages (which might dilute the effect), you still observe a statistically robust and consistent association.
# The effect persists across model specifications; The postSurge control confirms the model isn't confounded by extreme weather shocks
# The magnitude is modest but credible for environmental change processes. It is likley modest, because the impact of aquculture on salinity is more localised (in the absolute neighborhoods of the aquaculture ponds which is averaged across the village)


# 3. By State
regF3_1 <- feols(Saline_perc_norm  ~ Lag_Aqua * State | UniqueID + Year, 
                 data = data_filtered, vcov = ~UniqueID)
summary(regF3_1)

regF3_2 <- feols(Saline_perc_norm  ~ Lag_Aqua * State + postSurge | UniqueID + Year, 
                 data = data_filtered, vcov = ~UniqueID)
summary(regF3_2)

# Adding postSurge has no impact on the model or the coefficient of past aquaculture 

models <- list()
models[['Salinity (Village+Time FE)']] <- regF3_1
models[['Salinity (All controls+Village+Time FE)']] <- regF3_2
msummary(models, stars = c('*' = .1, '**' = .05, '***' = .01),gof_omit=c("BIC|AIC|RMSE|R2 Within Adj."),coef_omit=c("(Intercept)"), filename = 'table.rtf')

# Potential explanations for decreased salinity in Tn with increased aquaculture, but increased salinity in OD and AP
# Tamil Nadu has more freshwater resources than brackish water, but both play a significant role in its aquaculture industry. 
# While the state boasts 56,000 hectares of brackish water, it has a larger area of freshwater resources, including rivers, lakes, reservoirs, and tanks, totaling about 3.83 lakh hectares. 
# Brackish water aquaculture in Tamil Nadu, while increasing, is still at a smaller scale compared to freshwater aquaculture. 
# Source : https://tnfisheries.demodev.in/Aquaculture.html
# https://tnfisheries.demodev.in/includes/assets/cms_uploads/pdf/glance/FISHERIES_AT_A_GLANCE_2023-24_9604.pdf


# ----------------------
# VISUALIZATIONS
# ----------------------

# Basic Relationships 
a<-read.csv("visuals/AquaSalineSurge_Relationships.csv")

# Use standard errors to calculate CIs
a$upr <- a$est+(1.96*a$se)
a$lwr <- a$est-(1.96*a$se)

# Relevel the first grouping variable so that the plot is ordered the way I want
head(a)
a$Variable <- as.factor(a$Variable)
levels(a$Variable)
a$Variable <- factor(a$Variable, levels=c("Aquaculture (Past)","Salinity (Past)","Surge"))
a$Controls<-as.factor(a$Controls)

# Restrict the data only to estimates from select models (in this case all)
# a1 <- subset(a, Controls == "Village + Time FE" | Controls == "OD" | Controls == "AP" | Controls == "TN")
a1 <- subset(a, Controls == "Village + Time FE")


# Make sure that the Controls column is a factor and reorder as needed
levels(a1$Controls)
a1$Controls <- factor(a1$Controls, levels=c("Village + Time FE"))
a1$Controls<-as.factor(a1$Controls)


# Run the plot
# This is the plot with the horizontal lines
plot<-ggplot(a1,aes(est,Variable,group=Controls))+
  geom_errorbarh(aes(xmax = upr, xmin = lwr),height=0.2,position=position_dodge(width=0.2), colour="black")+
  geom_point(size=6, aes(shape=Controls,color=Controls),position=position_dodge(width=0.2))+
  theme_bw()+
  ggtitle("Estimated Effect on Aquaculture")+
  ylab("")+   theme(legend.position="right")+
  theme(legend.text = element_text(size = 12),legend.title = element_text(size = 14), legend.key.size = unit(1, 'cm'),legend.key.width = unit(0.4, "cm"),legend.spacing.x = unit(0.2, 'cm'))+
  theme(plot.title = element_text(hjust = 0),panel.border = element_rect(colour = "black",linewidth=1))+
  theme(axis.title = element_text(size = 12,face="bold"),axis.text=element_text(size=12,color="black"))+
  theme(plot.title = element_text(size = 14,face = "bold"))+xlab("") +
  theme(panel.grid.major = element_line(colour="grey", linewidth=0.5)) +
  theme(panel.grid.minor = element_line(colour="grey", linewidth=0.5))+
  theme(panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        panel.border = element_blank(),
        axis.line = element_line(color = 'black'))+theme(legend.title = element_blank(),
                                                         legend.key.size = unit(3, 'cm'))+
  geom_hline(yintercept = 4.5, linetype = "dotted",color ="black",linewidth = 1)+
  geom_hline(yintercept = 6.5, linetype = "dotted",color ="black",linewidth = 1)+
  geom_vline(xintercept = 0, linetype = "longdash",size = 1)+

  scale_color_manual(values = c( "#990000"))+
  scale_shape_manual(values = c(16,15,18,17))

print(plot)

ggsave("visuals/SetA_BasicRelationships.png",plot = plot, width = 8, height = 6, units = "in",dpi=600)



# Plot estimated effects for each state 
# Surge effect by state
# Extract coef and SE from regD1_2 for postSurge effects by State
coef_regD1_2 <- data.frame(
  State = c("Baseline (TN)", "AP", "OD"),
  Estimate = c(
    coef(regD1_2)["postSurge"],
    coef(regD1_2)["postSurge"] + coef(regD1_2)["postSurge:StateAP"],
    coef(regD1_2)["postSurge"] + coef(regD1_2)["postSurge:StateOD"]
  ),
  SE = c(
    sqrt(vcov(regD1_2)["postSurge", "postSurge"]),
    sqrt(vcov(regD1_2)["postSurge", "postSurge"] + 
           vcov(regD1_2)["postSurge:StateAP", "postSurge:StateAP"] +
           2 * vcov(regD1_2)["postSurge", "postSurge:StateAP"]),
    sqrt(vcov(regD1_2)["postSurge", "postSurge"] + 
           vcov(regD1_2)["postSurge:StateOD", "postSurge:StateOD"] +
           2 * vcov(regD1_2)["postSurge", "postSurge:StateOD"])
  )
)

coef_regD1_2 <- coef_regD1_2 %>%
  mutate(
    CI_lower = Estimate - 1.96 * SE,
    CI_upper = Estimate + 1.96 * SE
  )

ggplot(coef_regD1_2, aes(x = State, y = Estimate, ymin = CI_lower, ymax = CI_upper)) +
  geom_pointrange(color = "blue", size = 1.5) +
  geom_hline(yintercept = 0, linetype = "dashed", size = 1) +
  labs(
    title = "Effect of Surge on Aquaculture by States",
    y = "Estimated Effect (Coefficient)",
    x = "State"
  ) +
  theme_minimal(base_size = 16)

ggsave("visuals/SurgeEffect_ByState.png", width = 8, height = 6, units = "in",dpi=600)




# Salinity effect by State 
coef_regD1_4 <- data.frame(
  State = c("Baseline (TN)", "AP", "OD"),
  Estimate = c(
    coef(regD1_4)["avg_salinity_5yr"],
    coef(regD1_4)["avg_salinity_5yr"] + coef(regD1_4)["avg_salinity_5yr:StateAP"],
    coef(regD1_4)["avg_salinity_5yr"] + coef(regD1_4)["avg_salinity_5yr:StateOD"]
  ),
  SE = c(
    sqrt(vcov(regD1_4)["avg_salinity_5yr", "avg_salinity_5yr"]),
    sqrt(vcov(regD1_4)["avg_salinity_5yr", "avg_salinity_5yr"] + 
           vcov(regD1_4)["avg_salinity_5yr:StateAP", "avg_salinity_5yr:StateAP"] +
           2 * vcov(regD1_4)["avg_salinity_5yr", "avg_salinity_5yr:StateAP"]),
    sqrt(vcov(regD1_4)["avg_salinity_5yr", "avg_salinity_5yr"] + 
           vcov(regD1_4)["avg_salinity_5yr:StateOD", "avg_salinity_5yr:StateOD"] +
           2 * vcov(regD1_4)["avg_salinity_5yr", "avg_salinity_5yr:StateOD"])
  )
)

coef_regD1_4 <- coef_regD1_4 %>%
  mutate(
    CI_lower = Estimate - 1.96 * SE,
    CI_upper = Estimate + 1.96 * SE
  )

ggplot(coef_regD1_4, aes(x = State, y = Estimate, ymin = CI_lower, ymax = CI_upper)) +
  geom_pointrange(color = "orange", size = 1.5) +
  geom_hline(yintercept = 0, linetype = "dashed", size = 1) +
  labs(
    title = "Effect of Past Salinity on Aquaculture by States",
    y = "Estimated Effect (Coefficient)",
    x = "State"
  ) +
  theme_minimal(base_size = 16)

ggsave("visuals/SalinityEffect_ByState.png", width = 8, height = 6, units = "in",dpi=600)



# Feedback Salinity effect by state 

# Construct coefficient table for marginal effects of Lag_Aqua
coef_regF3_2 <- data.frame(
  State = c("Baseline (TN)", "AP", "OD"),
  Estimate = c(
    coef(regF3_2)["Lag_Aqua"],
    coef(regF3_2)["Lag_Aqua"] + coef(regF3_2)["Lag_Aqua:StateAP"],
    coef(regF3_2)["Lag_Aqua"] + coef(regF3_2)["Lag_Aqua:StateOD"]
  ),
  SE = c(
    sqrt(vcov(regF3_2)["Lag_Aqua", "Lag_Aqua"]),
    sqrt(
      vcov(regF3_2)["Lag_Aqua", "Lag_Aqua"] +
        vcov(regF3_2)["Lag_Aqua:StateAP", "Lag_Aqua:StateAP"] +
        2 * vcov(regF3_2)["Lag_Aqua", "Lag_Aqua:StateAP"]
    ),
    sqrt(
      vcov(regF3_2)["Lag_Aqua", "Lag_Aqua"] +
        vcov(regF3_2)["Lag_Aqua:StateOD", "Lag_Aqua:StateOD"] +
        2 * vcov(regF3_2)["Lag_Aqua", "Lag_Aqua:StateOD"]
    )
  )
)

coef_regF3_2 <- coef_regF3_2 %>%
  mutate(
    CI_lower = Estimate - 1.96 * SE,
    CI_upper = Estimate + 1.96 * SE
  )

# Plot
ggplot(coef_regF3_2, aes(x = State, y = Estimate, ymin = CI_lower, ymax = CI_upper)) +
  geom_pointrange(color = "firebrick", size = 1.5) +
  geom_hline(yintercept = 0, linetype = "dashed", size = 1) +
  labs(
    title = "Marginal effect of Past Aquaculture on Salinity by State",
    y = "Estimated Effect (Coefficient)",
    x = "State"
  ) +
  theme_minimal(base_size = 16)

ggsave("visuals/FeedbackSalinity_ByState.png", width = 8, height = 6, units = "in",dpi=600)




# 
# # ----------------------
# # OLDER REGRESSION MODELS
# # ----------------------
# 
# 
# # Model 1: Basic pooled OLS regression
# model_1 <- lm(Aqua_perc ~ postSurge + Saline_perc_norm + 
#                Sea_Dist + DEM_avg, data = data)
# summary(model_1)
# 
# # Model 2: Adding the interaction term
# model_2 <- lm(Aqua_perc ~ postSurge + Saline_perc_norm + 
#                (Saline_perc_norm*postSurge) + Sea_Dist + DEM_avg, data = data)
# summary(model_2)
# 
# 
# # Model 3 - State level FE + Two way clustering (also clustering for villages across years)
# # Relevel State so that TN is the reference category
# data$State <- relevel(factor(data$State), ref = "TN")
# 
# model_3 <- feols(Aqua_perc ~ postSurge + Saline_perc_norm + 
#                    (Saline_perc_norm*postSurge) + (State * postSurge) + DEM_avg | Sea_Dist, 
#                  data = data, 
#                  vcov = ~UniqueID ) 
# summary(model_3)
# 
# 
# models <- list()
# models[['Basic OLS']] <- model_1
# models[['OLS + Interaction']] <- model_2
# models[['Fixed Effects (Dist) + Clustered SE']] <- model_3
# msummary(models, stars = c('*' = .1, '**' = .05, '***' = .01),gof_omit=c("BIC|AIC|RMSE|R2 Within Adj."),coef_omit=c("(Intercept)"), filename = 'table.rtf')
# 
# # Interpretation: 
# # In TN, storm surges are associated with a +0.53 pp increase in aquaculture %, at average salinity and elevation.
# # A 1 SD increase in normalized salinity is associated with a +0.91 pp increase in aquaculture %.
# # The effect of salinity on aquaculture is less positive after a storm surge—by ~0.463 pp per SD.
# # In non-surge years, Andhra Pradesh has +0.66 pp higher aquaculture % than TN.
# # In non-surge years, Odisha has +0.21 pp higher aquaculture % than TN.
# # The storm surge effect in AP is 4.04 pp greater than in TN. Total effect = 0.53 + 4.04 ≈ +4.59 pp
# # The storm surge effect in OD is 1.04 pp less than in TN. Total effect = 0.53 − 1.04 ≈ −0.51 pp
# # Aquaculture slightly declines with elevation — by 0.06 pp for every 100 meters of elevation.
# # RMSE: ~5.2, meaning the average model prediction error is ~5 percentage points.
# # Adjusted R²: ~0.053 — low, but reasonable given village-level variation.
# # Within R²: ~0.05 — ~5% of within-village variation in aquaculture is explained.
# 
# # This suggests that storm surges are strongly associated with increased aquaculture in Andhra Pradesh, modest increases in Tamil Nadu, and even declines in Odisha. 
# # Higher salinity levels also generally promote aquaculture, but their positive impact weakens after a surge. 
# # Slight negative effects of elevation suggest topographic barriers to expansion. 
# # Overall, there are visible state-level differences in how aquaculture responds to climate shocks, shaped by both geography, policies, and prior land use dynamics.
# 
# 
# # For persistent Salinity - Using the 5-year average salinity
# # Model 4: Basic pooled OLS regression
# model_4 <- lm(Aqua_perc ~ postSurge + avg_salinity_5yr + 
#                 Sea_Dist + DEM_avg, data = data)
# summary(model_4)
# 
# # Model 5: Adding the interaction term
# model_5 <- lm(Aqua_perc ~ postSurge + avg_salinity_5yr + 
#                 avg_salinity_5yr*postSurge + Sea_Dist + 
#                 DEM_avg, data = data)
# summary(model_5)
# 
# 
# # Model 6: State level variations + clustering for villages across years
# model_6 <- feols(Aqua_perc ~ postSurge + avg_salinity_5yr + 
#                    (avg_salinity_5yr*postSurge) + (State * postSurge) + DEM_avg | Sea_Dist, 
#                  data = data, 
#                  vcov = ~UniqueID )
# summary(model_6)
# 
# models <- list()
# models[['Basic OLS']] <- model_4
# models[['OLS + Interaction']] <- model_5
# models[['Fixed Effects (Dist) + Clustered SE']] <- model_6
# msummary(models, stars = c('*' = .1, '**' = .05, '***' = .01),gof_omit=c("BIC|AIC|RMSE|R2 Within Adj."),coef_omit=c("(Intercept)"), filename = 'table.rtf')
# # Interpretation: 
# # Before a storm surge, a 1 unit increase in average salinity over the last 5 years is associated with a 1.44 percentage point increase in aquaculture.
# # suggests strong association of prior salinity levels with aquaculture.
# # After a storm surge, aquaculture increased by 0.503 pp on average in TN, holding average salinity constant.
# # StateAP: +0.665 => Andhra Pradesh villages have higher baseline aquaculture than TN, even before the storm.
# # StateOD: +0.218 => Odisha also shows higher aquaculture than TN pre-storm.
# # postSurge:StateAP: +4.04* => AP sees a much larger increase (4.04 pp) in aquaculture post-surge compared to TN.
# # postSurge:StateOD: –1.17* → Odisha sees a decline of ~1.17 pp relative to TN in post-surge aquaculture.
# # Interaction of salinity and storm : The effect of average salinity on aquaculture declines by 0.76 pp after a storm.
# # So post-surge, the marginal effect of salinity becomes: 1.44 – 0.76 ≈ 0.68 pp increase per unit salinity → salinity still increases aquaculture, but less strongly.
# 
# # Aquaculture is strongly predicted by long-term salinity, but that relationship weakens after storm surges. 
# # Andhra Pradesh shows much greater post-surge aquaculture expansion than Tamil Nadu, while Odisha shows a drop. 
# # Salinity, even averaged over 5 years, encourages aquaculture, but its effect is tempered after surge events, possibly due to damage, altered soil dynamics, or adaptation lag.
# 
# 
# 
# # With Salinity as binary 
# # Model 7: Basic pooled OLS regression
# model_7 <- lm(Aqua_perc ~ postSurge + Saline + 
#                 Sea_Dist + DEM_avg, data = data)
# summary(model_7)
# 
# # Model 8: Adding the interaction term
# model_8 <- lm(Aqua_perc ~ postSurge + Saline + 
#                 Saline * postSurge + Sea_Dist+ 
#                 DEM_avg, data = data)
# summary(model_8)
# 
# 
# # Model 9 - State level variations + FE for distance + clustering errors for villages across years
# model_9 <- feols(Aqua_perc ~ postSurge + Saline + 
#                    (Saline*postSurge) + (State * postSurge) + DEM_avg | Sea_Dist, 
#                  data = data, 
#                  vcov = ~UniqueID )
# summary(model_9)
# 
# models <- list()
# models[['Basic OLS']] <- model_7
# models[['OLS + Interaction']] <- model_8
# models[['Fixed Effects (State+Dist)']] <- model_9
# 
# msummary(models, stars = c('*' = .1, '**' = .05, '***' = .01),gof_omit=c("BIC|AIC|RMSE|R2 Within Adj."),coef_omit=c("(Intercept)"), filename = 'table.rtf')
# # Interpretation: 
# # In pre-surge periods, villages above the median salinity level have 0.6 percentage points more aquaculture than low-salinity villages.
# # PostSurge coefficient is for low-salinity village, i.e. In low-salinity villages, the storm surge is associated with a +0.20 pp increase in aquaculture. 
# # The lower statistical significance (p value 0.1%), implies not much evidence of change in low-saline villages after a surge.
# # Interaction between postSurge and Saline shows the additional effect of a storm surge in high-saline villages.
# # High-saline villages see an extra 0.52 pp increase in aquaculture after a storm compared to low-saline villages.
# # post-surge aquaculture in high-saline villages rises by: postSurge + postSurge:Saline = 0.23 + 0.52 = ~0.75 pp
# # Andhra Pradesh and Odisha both have higher baseline aquaculture (pre-surge) than Tamil Nadu (reference state), especially AP.
# # Huge post-surge increase in aquaculture in AP (vs TN): about +4.1 percentage points. (total = 4.1 + 0.23 = 4.33 pp)
# # Odisha sees a 0.92 pp smaller increase in aquaculture post-surge than TN. (total = 0.23 - 0.92 = -0.69 pp)
# # Higher elevation slightly reduces aquaculture by .07 pp per 100m elevation (likely due to lower access to brackish water).
# 
# # Before storms, aquaculture is concentrated in high-saline villages. 
# # After a storm surge, aquaculture rises more in these already saline areas, suggesting that saline villages may double down on aquaculture rather than revert to traditional agriculture. 
# # Low-saline villages show no significant response to storms. This reinforces the idea that salinity creates a path dependency for aquaculture expansion.
# # Andhra Pradesh continues to show the largest post-surge increases, hinting at policy support or infrastructure enabling rapid aquaculture expansion.
# # Odisha, by contrast, may face constraints or damage that depress post-surge growth. Could be owing to structural poverty. 
# 
# # Although, the r2 is less than it is with persistent salinity
# 
# 
# # Print results of different saline conditions 
# # All with - FE (Distance) + Clustered SEs + Interaction terms 
# models <- list()
# models[['Above median Saline']] <- model_9
# models[['Salinity Area']] <- model_3
# models[['Persistent Salinity']] <- model_6
# 
# lapply(models, function(m) names(coef(m)))
# 
# msummary(
#   models,
#   stars = c('*' = .1, '**' = .05, '***' = .01),
#   gof_omit = c("BIC|AIC|RMSE|R2 Within Adj."),
#   coef_omit = c("(Intercept)"),
#   coef_map = c(
#     "postSurge" = "Post-Surge (TN)",
#     "postSurge:StateAP" = "Post-surge (AP wrt TN)", 
#     "postSurge:StateOD" = "Post-surge (OD wrt TN)",
#     "Saline" = "Saline (above median)",
#     "Saline_perc_norm" = "Saline %",
#     "avg_salinity_5yr" = "Persistent Salinity (5-year avg)",
#     "postSurge:Saline" = "Post-Surge x Saline", 
#     "postSurge:Saline_perc_norm" = "Post-Surge x Saline %", 
#     "postSurge:avg_salinity_5yr" = "Post-Surge x Persistent Salinity", 
#     "StateAP" = "Pre-surge (AP wrt TN)",
#     "StateOD" = "Pre-surge (OD wrt TN)",
#     "DEM_avg" = "Elevation from the Sea"
#   ),
#   filename = 'table.rtf'
# )
# 
# 
# # Two-way fixed effects 
# # Note: Two-way fixed effect controls for village-level fixed effects (absorbing time-invariant differences between villages) - which is why NEAR_DIST and DEM_avg get dropped.
# # It controls for year fixed effects (accounting for changes common to all villages in a given year).
# # The model focuses on within-village variation over time rather than cross-sectional differences.
# 
# data_agg <- data %>%
#   group_by(UniqueID, Year) %>%
#   summarize(
#     Aqua_ha = mean(Aqua_ha, na.rm = TRUE),
#     Aqua_perc = mean(Aqua_perc, na.rm = TRUE),
#     Saline_perc_norm = mean(Saline_perc_norm, na.rm = TRUE),
#     postSurge = max(postSurge, na.rm = TRUE), # A village is affected if any record says so
#     salinity_storm_interaction <- Saline_perc_norm * postSurge,
#     DEM_avg = mean(DEM_avg, na.rm = TRUE),
#     NEAR_DIST = mean(NEAR_DIST, na.rm = TRUE),
# 
#     .groups = "drop"
#   )
# 
# summary(data_agg)
# 
# # Set up panel data structure
# pdata <- pdata.frame(data_agg, index = c("UniqueID", "Year"))
# 
# # Model 10: Fixed effects model (two-way fixed effect - village and year fixed effects)
# model_10 <- plm(Aqua_perc ~ postSurge + Saline_perc_norm + 
#                  postSurge * Saline_perc_norm, 
#                data = pdata, 
#                model = "within", 
#                effect = "twoways")
# summary(model_10)
# 
# # After a storm surge, aquaculture increases by ~0.28 pp, holding salinity constant.
# # In pre-surge years, a one-unit increase in normalized salinity reduces aquaculture by 0.27 pp.
# # After a storm, higher salinity is associated with more aquaculture: the negative effect of salinity flips post-surge.
# # => So, after a storm surge, salinity no longer deters aquaculture; in fact, the relationship flips to slightly positive, though small.
# # After storm-induced salinization, farmers may switch to aquaculture in saline lands.
# 
# # The R2 may seem low, but it is expected because two-way FE models absorb a lot of variation through FE. 
# # Besides, this model explain within-village, within-year variation, and not between village or across year variations. 
# # Patterns in coefficients being significant is key here. 
# 
# # After a storm surge, salinity becomes a trigger for aquaculture — possibly because traditional crops can no longer grow. 
# # This supports the idea that post-surge salinization causes land use change toward aquaculture, particularly in areas that become newly saline.
# 
# 
# # Testing for serial correlation
# pbgtest(model_10)
# 
# # This tests whether the error terms within a unit over time (e.g., for a village) are correlated.
# # Serial correlation violates a key assumption of standard regression models
# 
# # Since serial correlation is present, we could: 
# # 1. Use robust standard errors with vcovHC()
# # 2. Consider including lagged dependent variables in the model 
# # 3. Consider using first differences or other dynamic panel data methods
# # 4. Try Newey-West standard errors or other corrections
# 
# # Using robust standard errors
# robust_model_10 <- coeftest(model_10, vcov = vcovHC(model_10, type = "HC1"))
# print(robust_model_10)
# # The effect estimates remain mostly the same and highly statistically significant, despite accounting for serial correlation and heteroskedasticity. 
# 
# 
# models <- list()
# models[['Robust Two-way FE']] <- robust_model_10
# 
# msummary(models, stars = c('*' = .1, '**' = .05, '***' = .01),gof_omit=c("BIC|AIC|RMSE|R2 Within Adj."),coef_omit=c("(Intercept)"), filename = 'table.rtf')
# # All three effects remain statistically significant at the 0.001 level, even with robust standard errors.
# # Before surge: Salinity reduces aquaculture — each unit of Saline_perc_norm is associated with ~0.27 pp lower aquaculture.
# # After surge: The interaction offsets this effect — total effect of salinity becomes:
# # − 0.266 +  0.275=  +  0.009
# # No negative effect anymore — salinity no longer deters aquaculture.
# 
# # => Villages that become saline after a storm are more likely to switch to aquaculture.
# # The interaction effect is both statistically and substantively important, even after accounting for village and year fixed effects and correcting for serial correlation.
# 
# 
# 
# 
# 
# 
# # ----------------------
# # VISUALIZATIONS
# # ----------------------
# 
# # Model 3 (Saline Area %) visualisation - moving FE to a predictor to enable prediction 
# data$State <- relevel(factor(data$State), ref = "TN")
# model_3 <- feols(Aqua_perc ~ postSurge + Saline_perc_norm + 
#                    (Saline_perc_norm*postSurge) + (State * postSurge) + DEM_avg + Sea_Dist, 
#                  data = data, 
#                  vcov = ~UniqueID ) 
# 
# # Visualizing Predicted Values by State and postSurge
# # Step 1: Create a new data frame for prediction
# # Get means of continuous covariates
# mean_salinity <- mean(data$Saline_perc_norm, na.rm = TRUE)
# mean_dem <- mean(data$DEM_avg, na.rm = TRUE)
# 
# # Create prediction grid
# pred_data <- expand.grid(
#   postSurge = c(0, 1),
#   State = c("AP", "OD", "TN"),
#   Saline_perc_norm = mean_salinity,
#   DEM_avg = mean_dem,
#   Sea_Dist = unique(data$Sea_Dist)[1]  # pick one for fixed effects
#   )
# 
# # Step 2: Generate predicted values using predict()
# # Generate predicted values (with confidence intervals)
# pred_data$pred <- predict(model_3, newdata = pred_data, se = TRUE)$fit
# pred_data$se <- predict(model_3, newdata = pred_data, se = TRUE)$se.fit
# 
# # Confidence intervals
# pred_data$lower <- pred_data$pred - 1.96 * pred_data$se
# pred_data$upper <- pred_data$pred + 1.96 * pred_data$se
# 
# 
# # Step 3: Plot the predicted values
# ggplot(pred_data, aes(x = factor(postSurge), y = pred, fill = State)) +
#   geom_col(position = position_dodge(0.6), width = 0.5) +
#   geom_errorbar(aes(ymin = lower, ymax = upper), position = position_dodge(0.6), width = 0.2) +
#   labs(
#     x = "Post Surge (0 = Before, 1 = After)",
#     y = "Predicted Aquaculture %",
#     title = "Predicted Aquaculture Share by State and Storm Surge",
#     fill = "State"
#   ) +
#   theme_minimal()
# 
# 
# 
# # Robust Model 10 Visualization 
# # Plot predicted Aqua_perc across a range of Saline_perc_norm values using two-way FE (village and year) model. 
# # In regression analysis, especially when using models with interactions or nonlinearities, the marginal effect tells us: 
# # "How much does the dependent variable change when one independent variable increases by a small amount, holding everything else constant?"
# # In this case, since we are using an interaction in the equation: 
# # Aqua_perc = b0 + b1 Saline_perc_norm + b2 postSurge + b3 (postSurge * Saline_perc_norm) 
# # dAqua_perc/dsaline = b1 + b3 x postSurge (differential eqn)
# # then b1 is the baseline sensitivity of aquaculture to salinity and b3 is additional sensitivity due to storm 
# # Visualizing this : 
# 
# # Step 1: Generate a smooth sequence of salinity values and create prediction data by State and Surge status
# sal_seq <- seq(min(data$Saline_perc_norm, na.rm = TRUE),
#                max(data$Saline_perc_norm, na.rm = TRUE), 
#                length.out = 100)
# 
# # Create all combinations of salinity x postSurge x State
# states <- unique(data$State)
# marg_data <- expand.grid(
#   Saline_perc_norm = sal_seq,
#   postSurge = c(0, 1),
#   State = states
# )
# 
# # Add interaction term
# marg_data$postSurge_Saline <- marg_data$postSurge * marg_data$Saline_perc_norm
# 
# # Step 3: Get robust coefficients and variance-covariance matrix
# coefs <- coef(robust_model_10) |> as.numeric()
# vcov_mat <- vcovHC(model_10, type = "HC1")
# 
# 
# # Design matrix
# X <- model.matrix(~ postSurge + Saline_perc_norm + postSurge_Saline - 1, data = marg_data)
# # The -1 removes the intercept
# 
# dim(X)        # (n, k)
# length(coefs) # should be k
# 
# 
# # Predictions and standard errors
# marg_data$predicted <- as.numeric(X %*% coefs)
# marg_data$se <- sqrt(rowSums((X %*% vcov_mat) * X))
# 
# 
# # Confidence intervals (95%)
# marg_data <- marg_data %>%
#   mutate(
#     lower = predicted - 1.96 * se,
#     upper = predicted + 1.96 * se,
#     SurgePeriod = ifelse(postSurge == 1, "Post-surge", "Pre-surge")
#   )
# 
# # Step 4: Plot
# ggplot(marg_data, aes(x = Saline_perc_norm, y = predicted, color = SurgePeriod, fill = SurgePeriod)) +
#   geom_line(size = 1.2) +
#   geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.2, linetype = 0) +
#   facet_wrap(~ State) +
#   labs(
#     title = "Marginal Effect of Salinity on Aquaculture by State",
#     x = "Normalized Salinity (%)",
#     y = "Predicted Aquaculture (%)",
#     color = "Storm Surge Period",
#     fill = "Storm Surge Period"
#   ) +
#   theme_minimal(base_size = 14)
# 
# # The downward pre-surge suggests negative marginal effect of salinity 
# # the upward post-surge suggests that surges reverse the penalty of salinity on aquaculture 
# # The gap between lines increases with salinity → storm surges shift the relationship especially for more saline villages.
# 
# 
# # Option 2: plot group means instead of fitted lines 
# # Bin salinity for grouping (e.g., 50 bins across full range)
# data_binned <- data %>%
#   mutate(
#     Salinity_bin = cut(Saline_perc_norm, breaks = 10, include.lowest = TRUE),
#     SurgePeriod = ifelse(postSurge == 1, "Post-surge", "Pre-surge")
#   ) %>%
#   group_by(Salinity_bin, postSurge, SurgePeriod, State) %>%
#   summarise(
#     Saline_perc_norm = mean(Saline_perc_norm, na.rm = TRUE),  # bin center
#     mean_aqua = mean(Aqua_perc, na.rm = TRUE),
#     se = sd(Aqua_perc, na.rm = TRUE) / sqrt(n()),
#     n = n(),
#     .groups = "drop"
#   ) %>%
#   mutate(
#     lower = mean_aqua - 1.96 * se,
#     upper = mean_aqua + 1.96 * se
#   )
# 
# # Calculate mean salinity (for vertical line)
# sal_mean <- mean(data$Saline_perc_norm, na.rm = TRUE)
# 
# # Plot group means with CI ribbons
# ggplot(data_binned, aes(x = Saline_perc_norm, y = mean_aqua,
#                         color = SurgePeriod, fill = SurgePeriod)) +
#   geom_line(size = 1.2) +
#   geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.2, linetype = 0) +
#   geom_vline(xintercept = sal_mean, linetype = "dashed", color = "black") +
#   facet_wrap(~ State) +
#   labs(
#     title = "Predicted Aquaculture by Salinity, State and Storm Surge",
#     x = "Normalized Salinity (%)",
#     y = "Predicted Aquaculture (%)",
#     color = "Storm Surge Period",
#     fill = "Storm Surge Period"
#   ) +
#   theme_minimal(base_size = 14)
# 
# 
# 
# 
# 
# 
# # Compare the first set of regression coefficients 
# a<-read.csv("outputs/RQ1_2013-25_AquaHa_Storm_Salinity_Models.csv")
# 
# # Use standard errors to calculate CIs
# a$upr <- a$est+(1.96*a$se)
# a$lwr <- a$est-(1.96*a$se)
# 
# # Relevel the first grouping variable so that the plot is ordered the way I want
# head(a)
# a$Variable <- as.factor(a$Variable)
# levels(a$Variable)
# a$Variable <- factor(a$Variable, levels=c("Salinity (5 yrs)","Saline Area", "Surge"))
# a$Controls<-as.factor(a$Controls)
# 
# 
# # Restrict the data only to estimates from select models
# a1 <- subset(a, Controls == "OLS + Interaction" | Controls == "FE (State + Sea Distance)" | Controls == "Persistent Salinity (OLS + Interaction)" | Controls == "Persistent Salinity (FE)")
# 
# # Make sure that the Controls column is a factor and reorder as needed
# levels(a1$Controls)
# a1$Controls <- factor(a1$Controls, levels=c("OLS + Interaction","FE (State + Sea Distance)","Persistent Salinity (OLS + Interaction)", "Persistent Salinity (FE)"))
# a1$Controls<-as.factor(a1$Controls)
# 
# # Run the plot
# # This is the plot with the horizontal lines
# plot<-ggplot(a1,aes(est,Variable,group=Controls))+
#   geom_errorbarh(aes(xmax = upr, xmin = lwr),height=0.5,position=position_dodge(width=0.5), colour="black")+
#   geom_point(size=5, aes(shape=Controls,color=Controls),position=position_dodge(width=0.5))+
#   theme_bw()+
#   ggtitle("Relationship between Aquaculture, Surge & Salinity across model specs")+
#   ylab("")+   theme(legend.position="right")+
#   theme(legend.text = element_text(size = 12),legend.title = element_text(size = 12))+
#   theme(plot.title = element_text(hjust = 0.5),panel.border = element_rect(colour = "black",linewidth=1))+
#   theme(axis.title = element_text(size = 16,face="bold"),axis.text=element_text(size=16,color="black"))+
#   theme(plot.title = element_text(size = 16,face = "bold"))+xlab("") +
#   theme(panel.grid.major = element_line(colour="grey", linewidth=0.5)) +
#   theme(panel.grid.minor = element_line(colour="grey", linewidth=0.5))+
#   theme(panel.grid.minor = element_blank(),
#         panel.grid.major = element_blank(),
#         panel.border = element_blank(),
#         axis.line = element_line(color = 'black'))+theme(legend.title = element_blank(), 
#                                                          legend.key.size = unit(3, 'cm'))+
#   geom_hline(yintercept = 4.5, linetype = "dotted",color ="black",linewidth = 1)+
#   geom_hline(yintercept = 6.5, linetype = "dotted",color ="black",linewidth = 1)+
#   geom_vline(xintercept = 0, linetype = "longdash",size = 0.5)+
#   
#   scale_color_manual(values = c( "#ffab33","#990000", "#666633", "#003300"))+
#   scale_shape_manual(values = c(16,15,18,17))
# 
# plot(plot)
# 
# ggsave("outputs/RQ1_2013-25_AquaHa_Storm_Salinity_Models.png",plot = plot, width = 12, height = 6, units = "in",dpi=300)
# 
# 
# 
# 
# 
# # Visualize the interaction effect
# # Create predicted values based on Model 2
# grid_data <- expand.grid(
#   Saline_perc_norm = seq(min(data$Saline_perc_norm, na.rm = TRUE),
#                             max(data$Saline_perc_norm, na.rm = TRUE),
#                             length.out = 100),
#   postSurge = c(0, 1),
#   NEAR_DIST = mean(data$NEAR_DIST, na.rm = TRUE),
#   DEM_avg= mean(data$DEM_avg, na.rm = TRUE)
# )
# 
# grid_data$salinity_storm_interaction <- grid_data$Saline_perc_norm * grid_data$postSurge
# grid_data$predicted <- predict(model_2, newdata = grid_data)
# 
# # Plot the interaction effect
# ggplot(grid_data, aes(x = Saline_perc_norm, y = predicted, color = factor(postSurge))) +
#   geom_line() +
#   labs(title = "Interaction Effect of Salinity and Storm Surge on Aquaculture",
#        x = "Saline Area Percentage",
#        y = "Predicted Aquaculture",
#        color = "Storm Surge") +
#   theme_minimal()
# 
# 
# 
# # Map Flood - Saline areas 
# # map1 <- merge(map_dist_shp, dist_state_data, by="District_name", all.x=T)
# # 
# # # Poverty 
# # summary(map1$poverty_share) 
# # # Some districts are outliers i.e. have too high poverty (>0.8) may want to remove those for mapping? 
# # map_poverty <- ggplot() +
# #   geom_sf(data=map_dist_shp, colour = "lightgrey", fill="lightgrey")+ # To map all districts even if data not available
# #   geom_sf(data=map1, aes(fill = poverty_share), colour=alpha("dark grey", 1/2), size=0.2) +
# #   geom_sf(data=map_state_shp, colour="white", fill=NA, size = 4) + # draw polygon outline of states and fill nothing
# #   coord_sf() +
# #   theme_minimal() +
# #   ggtitle("Poverty share across districts: % of people in 0-20% MPCE Quantile") +
# #   theme(axis.line = element_blank(), axis.text = element_blank(),
# #         axis.ticks = element_blank(), axis.title = element_blank(), 
# #         plot.title = element_text(size = 12),
# #         legend.key.size = unit(2, "lines"),
# #         panel.grid.major = element_blank(),
# #         panel.grid.minor = element_blank()) +
# #   scale_fill_viridis(option="cividis", breaks = seq(0, 1, by = 0.2))
# # # Adjust legend position and size
# # map_poverty <- map_poverty +
# #   theme(legend.position = "bottom",  # Adjust legend position
# #         legend.text = element_text(size = 10),  # Adjust legend text size
# #         legend.title = element_text(size = 10))  # Adjust legend title size
# # #Export map
# # ggsave("Map_poverty.png",plot = map_poverty, width = 8, height = 8, units = "in",dpi=600)
# # 
