# RQ1: In what ways is the aquaculture land transition related to salinity prevalence and history of storm surges,
# accounting for place-invariant factors such as distance/elevation from the sea or soil conditions, or time based shocks such as policy changes, inflation, or market conditions

# RQ2: How do these relationships vary by states? 

# Optional RQ3: Is there a predictive lag between surge impact and land change outcome? 


# Aquaculture Analysis: Panel Data Regression Models
# This script analyzes the relationship between aquaculture land change, storm surges, salinity, and other socio-environmental factors

# As Moran's I is positive, we need to also account for spatial autocorrelation, ie villages near each other tend to have simlar levels of aquaculture => the residuals in the regression may also be spatially autocorrelated violating the classical assumption about independence. 
# vcov at the village level corrects for within-unit correlation over time (serial correlation), and heteroskedasticity within village clusters 
# it does not adjust for spatial dependence between nearby villages
# Try spatially clustered standard errors - Conley Standard Errors that adjust for spatial correlation across observations within a specified distance (5km, 10km...)


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
library(slider)
library(margins) # For visualizing 
library(sf)
library(spdep) #for conley spatial standard error 
library(spatialreg) #for conley spatial standard error 

data <- read.csv(unz("data/aqua_salinity_surge_1990-2025.zip", 
                           "aqua_salinity_surge_1990-2025.csv"))

# Need lat long in degrees for Conley
village_geometry <- st_read("data/shp/States_AllData_geometry.shp")
village_wgs <- st_transform(village_geometry, crs = 4326)

# Check for invalid geometries
sum(!st_is_valid(village_wgs))  # number of invalid geometries

village_wgs <- village_wgs %>%
  mutate(geometry = st_make_valid(geometry))

village_wgs <- village_wgs %>%
  mutate(
    Longitude_d = st_coordinates(st_centroid(geometry))[, 1],
    Latitude_d = st_coordinates(st_centroid(geometry))[, 2]
  )

head(village_wgs)

degrees <- village_wgs %>%
  select(UniqueID, Longitude_d, Latitude_d)

data <- data %>%
  left_join(degrees, by = "UniqueID")
  
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

head(data)

# Add a flag variable for any key missing variables 
# If any of these is NA, make flag = NA : Aqua_perc , postSurge , avg_salinity_5yr , Lag_Aqua
data$flag <- ifelse(rowSums(is.na(data[c("Aqua_perc", "postSurge", "avg_salinity_5yr", "Lag_Aqua")])) > 0, NA, 0)

table(data$flag, useNA = "always")  # Shows count of flagged vs non-flagged rows
sum(is.na(data$flag))   # Total number of rows with missing variables


# ---------------------------------------
# REGRESSION MODELS : BASIC RELATIONSHIPS
# ---------------------------------------
# Two-way fixed effect models 

# SET A : Without and with controls 
# 1. Y = Salinity, X = surge (without controls except time and unit)
regA1_1 <- feols(avg_salinity_5yr ~ postSurge + flag | Year, 
                data = data, vcov = ~UniqueID)
summary(regA1_1)

regA1_2 <- feols(avg_salinity_5yr ~ postSurge + flag | UniqueID + Year, 
                data = data, vcov = ~UniqueID)
summary(regA1_2)

# Time FE: 
# Holding the year constant, when comparing across places, Villages in post-surge years have, on average, 0.608 units higher 5-year average salinity, compared to not affected places. 
# Note that this is a pooled estimate across all villages without controlling for baseline salinity differences between villages.
# This only explains 10% of the variation, and RMSE = 0.8 implies prediction errors are fairly large. 

# Village + Time FE
# This controls for any time-invariant village-level characteristics, like soil type, geographical location with respect to the sea, persistent farming practices, etc.
# When comparing villages to themselves overtime, the relationship of salinity with surges reverses, that is within a given village, salinity is lower by 0.066 units in post-surge periods compared to pre-surge periods.
# Within R² = 0.00022: only a tiny portion of within-village variation in salinity is explained by postSurge.
# This suggests that there maybe systematic differences between places that are affected vs not affected by storms. Those affected by storms are likley to have a higher salinity baseline. 
# while places that are affected by surges have a higher baseline salinity (i.e. are more likely to be saline), they experience relatively less salinity after being affected by storms. 

# The flip suggests that areas that get surges more often may already have higher salinity. The negative effect may indicate a temporary leaching/flushing effect or behavioral change (e.g., fallowing land post-surge) or higher soil moisture

# apply spatial standard error 

regA1_2c <- feols(avg_salinity_5yr ~ postSurge + flag | UniqueID + Year, 
                 data = data)
summary(regA1_2c)

se_conley <- vcov_conley(
  regA1_2c,
  lat = ~Latitude_d,
  lon = ~Longitude_d,
  cutoff = 5  # in kilometers
)
etable(regA1_2c, vcov = se_conley)

# The coefficient is not statistically significant, ie, the ratio of the coefficient to its SE (t-statistic) is less than expected (1.65 for 10%, 1.96 for 5%, and 2.58 for 1% significance) 
# This is mainly because the effect is too small wrt to the SE 
# Therefore, we cannot reject the null hypothesis that postSurge has no effect on avg_salinity_5yr.
# See geographical heterogeneity in this - AP is significant at 1%, TN at 5% and OD not significantly different from TN. 



# 2. Y = Aquaculture, X = Surge (without controls except time and unit, with controls)
regA2_1 <- feols(Aqua_perc ~ postSurge + flag | Year, 
                 data = data, vcov = ~UniqueID)
summary(regA2_1)

regA2_2 <- feols(Aqua_perc ~ postSurge + flag | UniqueID + Year, 
                 data = data, vcov = ~UniqueID)
summary(regA2_2)

# Time FE: 
# Across all villages (this includes both within- and between-village variation), in years after a surge, the average aquaculture share is nearly 1 percentage point higher, compared to pre-surge years.
# Time + Villge FE: 
# Within a village over time, aquaculture increases by ~0.23 percentage points after a surge, holding constant all time-invariant village characteristics and year effects.
# This is accounting for place-based time-invariant village characteristics and time shocks that affect all villages (policies, inflation, prices, etc.), and therefore likely a causal effect.  

# After accounting for time-invariant differences across villages, storm surges are still positively associated with more aquaculture, but the effect size drops a lot — from 0.92% to 0.23%.
# This pattern is consistent with adaptation behavior: Surge hits → People shift toward aquaculture
# But the very low within R² suggests surge is only one small part of the story — most of the variation in aquaculture adoption is driven by other factors (e.g., market access, policy, infrastructure, long-run salinity trends, landholding).


# apply spatial standard error 
regA2_2c <- feols(Aqua_perc ~ postSurge + flag | UniqueID + Year, 
                 data = data)
summary(regA2_2c)

se_conley <- vcov_conley(
  regA2_2c,
  lat = ~Latitude_d,
  lon = ~Longitude_d,
  cutoff = 5  # in kilometers
)
etable(regA2_2c, vcov = se_conley)

# Although the effect size is large (23.1%) the ratio with SE (0.1887) is ~1.224 which is well below the 1.645 for 10% significance 
# Cannot exclude the null hypothesis that surge has no effect on aquaculture uptake 


# 3. Y = Aquaculture, X = Salinity in previous period (without controls except time and unit, with controls)
regA3_1 <- feols(Aqua_perc ~ avg_salinity_5yr + flag | Year, 
                 data = data, vcov = ~UniqueID)
summary(regA3_1)

regA3_2 <- feols(Aqua_perc ~ avg_salinity_5yr + flag | UniqueID + Year, 
                 data = data, vcov = ~UniqueID)
summary(regA3_2)

# Time FE: 
# 1-unit increase in 5-year salinity is associated with a 0.96 percentage point increase in aquaculture share across villages.
# higher salinity is correlated with more aquaculture at the aggregate level, but doesn't control for village-specific unobservables.

# Time and village FE: 
# 1-unit increase in salinity over the past 5 years is associated with a 0.2 percentage point increase in aquaculture share within the same village.
# Increasing salinity leads to relatively small but significant aquaculture expansion over time, within villages.

# These models provide evidence for adaptive behavior: Villages expand aquaculture in response to worsening salinity.
# The fact that the association holds even after controlling for fixed village characteristics, strengthens the interpretation that salinity increases precede aquaculture growth
# Although the effect size is modest => effect size is modest, and most of the change in aquaculture is driven by other time-varying factors.


# apply spatial standard error 
regA3_2c <- feols(Aqua_perc ~ avg_salinity_5yr + flag | UniqueID + Year, 
                 data = data)
summary(regA3_2c)

se_conley <- vcov_conley(
  regA3_2c,
  lat = ~Latitude_d,
  lon = ~Longitude_d,
  cutoff = 5  # in kilometers
)
etable(regA3_2c, vcov = se_conley)

# # Coefficient = 0.1988 - >  1-unit increase in average salinity over 5 years is associated with a 20.34 percentage point increase in aquaculture share (Aqua_perc).
# This suggests a strong positive association between higher long-term salinity and aquaculture area.
# # Standard Error = 0.1035 -> With a t-statistic of ~1.92 (0.1988 / 0.1035), this result is statistically significant at the 10% level, denoted by ..



# 4. Y = Aquaculture, X = Aqua time - 1 (without controls except time and unit, with controls)
# to assess persistence or path dependence in aquaculture land use decisions — whether villages that adopt aquaculture continue doing so.

regA4_1 <- feols(Aqua_perc ~ Lag_Aqua + flag | Year, 
                 data = data, vcov = ~UniqueID)
summary(regA4_1)

regA4_2 <- feols(Aqua_perc ~ Lag_Aqua + flag | UniqueID + Year, 
                 data = data, vcov = ~UniqueID)
summary(regA4_2)

# Aquaculture in the previous period remains consistently high and significant predictor or aquaculture in the current period

# Time FE: 
# One percentage point increase in aquaculture in year t–1 is associated with a 0.87 percentage point increase in year t. This implies very strong persistence of aquaculture from year to year 
# Nearly all the explanatory power comes from within-village year-to-year changes.

# Time and Village FE: 
# within the same village, aquaculture in year t is strongly predicted by aquaculture in year t–1 — a 0.64 percentage point increase per 1 unit aquaculture in the previous year.
# Within R² = 0.38: even after accounting for all fixed effects, the lag explains 38% of within-village variation, which is substantial.

# This suggests path dependence in aquaculture land use even after accounting for permanent village traits (e.g., location, infrastructure, soil).
# Aquaculture is highly persistent — villages that adopted it in the past tend to maintain or expand it.
# This persistence holds even when controlling for time-invariant village characteristics, indicating path dependency and potential lock-in once aquaculture is adopted.
# This may reflect: Economic investments (e.g., ponds, labor skills), Institutional or policy factors, Environmental feedbacks (e.g., increased salinity reinforcing suitability)

# apply spatial standard error 
regA4_2c <- feols(Aqua_perc ~ Lag_Aqua + flag | UniqueID + Year, 
                 data = data)
summary(regA4_2c)

se_conley <- vcov_conley(
  regA4_2c,
  lat = ~Latitude_d,
  lon = ~Longitude_d,
  cutoff = 5  # in kilometers
)
etable(regA4_2c, vcov = se_conley)

# Standard Error = 0.015, p < 0.001 -> The estimate is highly statistically significant. 
# strong evidence of path dependence or inertia in aquaculture. Once adopted, aquaculture practices tend to persist, possibly due to:
# Sunk investments in ponds or infrastructure, Skill specialization or market commitments, Environmental lock-ins (e.g., salinized land).



# 5. Y = Aquaculture, X = Surge, X = Salinity, X = Aqua-1 (without controls except time and unit, with controls)
regA5_1 <- feols(Aqua_perc ~ postSurge + avg_salinity_5yr + Lag_Aqua + flag | Year, 
                 data = data, vcov = ~UniqueID)
summary(regA5_1)

regA5_2 <- feols(Aqua_perc ~ postSurge + avg_salinity_5yr + Lag_Aqua + flag | UniqueID + Year, 
                 data = data, vcov = ~UniqueID)
summary(regA5_2)

# When including all predictors, they all seem to have a strong and significant relationship with aquaculture, both within and across comparisons. 
# Time FE: 
# Past aquaculture has strong persistence effect: Strong persistence — a 1% point increase in aquaculture last year is associated with 0.87% point increase this year.
# Past salinity is associated with an increase in aquaculture: +0.15% point increase in aquaculture.
# Surge years are associated with 0.09% point higher aquaculture, on average.

# Time and Village FE: 
# Persistence of aquaculture remains strong even after accounting for village-level confounders - 1 unit in past aquaculture accounts for 0.638 unit increase in future aquaculture 
# Surge is associated with a +0.12% point increase in aquaculture.
# Salinity increases are positively associated with +0.07% point increase in aquaculture.
# These results provide robust evidence for an adaptive feedback loop: Storm surge => Salinity increases => Aquaculture increases => Persistence over time

# apply spatial standard error 
regA5_2c <- feols(Aqua_perc ~ postSurge + avg_salinity_5yr + Lag_Aqua + flag | UniqueID + Year, 
                 data = data)
summary(regA5_2c)

se_conley <- vcov_conley(
  regA5_2c,
  lat = ~Latitude_d,
  lon = ~Longitude_d,
  cutoff = 5  # in kilometers
)
etable(regA5_2c, vcov = se_conley)

# Slight increase in aquaculture after surge, weakly significant.
# Higher salinity over past 5 years is weakly associated with more aquaculture.
# Very strong and highly significant autocorrelation — villages with aquaculture in t−1 tend to continue it in t.


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







# Set B : With salinity and surge interactions 
# To assess if there is a moderating effect of surge on the salinity and aquaculture relationship.
# 1. Y = Aquaculture, X = Surge, X = Salinity, X = Aqua-1 

regB1_1 <- feols(Aqua_perc ~ postSurge + avg_salinity_5yr + Lag_Aqua + postSurge * avg_salinity_5yr + flag | Year, 
                 data = data, vcov = ~UniqueID)
summary(regB1_1)

regB1_2 <- feols(Aqua_perc ~ postSurge + avg_salinity_5yr + Lag_Aqua + postSurge * avg_salinity_5yr + flag | UniqueID + Year, 
                 data = data, vcov = ~UniqueID)
summary(regB1_2)

# Time FE
# In non-saline areas, a surge year increases aquaculture by ~0.11% point.
# In non-surge years, a 1-unit increase in salinity increases aquaculture by ~0.24% point.
# After a surge, the effect of salinity on aquaculture is lower by 0.16 → Total effect of salinity = 0.24 – 0.16 = ~0.08.
# Salinity is positively associated with aquaculture in normal years.
# After a surge, the strength of that association weakens — perhaps because surge acts as a shock to the existing aquaculture practices or immediate recovery efforts shift behavior.

# Time and village FE: 
# Surges increase aquaculture in low-salinity areas - In villages with zero salinity, surges are associated with a +0.12% point increase in aquaculture.
# Salinity alone does not predict aquaculture in normal years (not significant) - In non-surge years, salinity has a small and statistically insignificant effect on aquaculture.
# Surges amplify salinity-driven aquaculture expansion - After a surge, the effect of salinity increases by ~0.09% point per unit of salinity.
# In normal years, salinity does not have an impact on aquaculture uptake, however, after a surge, salinity becomes a stronger driver of aquaculture.
# This points to interactive adaptation: when a surge occurs in already saline areas, people are more likely to expand aquaculture.
# In villages already facing salinity, storm surges amplify the move toward aquaculture.
# It is not just the surge or salinity alone, but salinity post-surge that accelerates aquaculture expansion.

# apply spatial standard error 
regB1_2c <- feols(Aqua_perc ~ postSurge + avg_salinity_5yr + Lag_Aqua + postSurge * avg_salinity_5yr + flag | UniqueID + Year, 
                 data = data)
summary(regB1_2c)

se_conley <- vcov_conley(
  regB1_2c,
  lat = ~Latitude_d,
  lon = ~Longitude_d,
  cutoff = 5  # in kilometers
)
etable(regB1_2c, vcov = se_conley)

# Only lag aqua is significant. But also, within R2 suggests that 38% of within village variation is captured by these three measures.  




models <- list()
models[['AC (All controls+Time FE)']] <- regB1_1
models[['AC (All controls+Village+Time FE)']] <- regB1_2

msummary(models, stars = c('*' = .1, '**' = .05, '***' = .01),gof_omit=c("BIC|AIC|RMSE|R2 Within Adj."),coef_omit=c("(Intercept)"), filename = 'table.rtf')



models <- list()
models[['Salinity (1)']] <- regA1_2

models[['Aquaculture (1)']] <- regA3_2

models[['Aquaculture (2)']] <- regA2_2

models[['Aquaculture (3)']] <- regA4_2

models[['Aquaculture (4-All)']] <- regA5_2

models[['Aquaculture (5-Interaction)']] <- regB1_2

msummary(models, stars = c('*' = .1, '**' = .05, '***' = .01),gof_omit=c("BIC|AIC|RMSE|R2 Within Adj."),coef_omit=c("(Intercept)"), filename = 'table.rtf')




# Set C : State variation

# 1. Y = Salinity, X = surge (without controls except time and unit, with controls)
regC1_1 <- feols(avg_salinity_5yr ~ postSurge + postSurge * State + flag | Year, 
                 data = data, vcov = ~UniqueID)
summary(regC1_1)

regC1_2 <- feols(avg_salinity_5yr ~ postSurge + postSurge * State + flag | UniqueID + Year, 
                 data = data, vcov = ~UniqueID)
summary(regC1_2)



# apply spatial standard error 
regC1_2c <- feols(avg_salinity_5yr ~ postSurge + postSurge * State + flag | UniqueID + Year, 
                 data = data)
summary(regC1_2c)

se_conley <- vcov_conley(
  regC1_2c,
  lat = ~Latitude_d,
  lon = ~Longitude_d,
  cutoff = 5  # in kilometers
)
etable(regC1_2c, vcov = se_conley)

# After surge events, average salinity decreases by about 0.21 units in TN. Significant at 5% => In Tamil Nadu, post-surge years are associated with a significant decrease in 5-year salinity
# In Andhra Pradesh (AP), the effect of surge on salinity is 0.4784 units higher than Tamil Nadu (0.4784-0.2102 = +0.2682,). Highly significant at 1% 
# In Odisha (OD), the surge effect on salinity is not significantly different from Tamil Nadu. => In OD, the effect is +0.113 higher than TN, i.e., net effect = -0.0973, not statistically significant



# 2. Y = Aquaculture, X = Surge (without controls except time and unit, with controls)
regC2_1 <- feols(Aqua_perc ~ postSurge + postSurge * State + flag | Year, 
                 data = data, vcov = ~UniqueID)
summary(regC2_1)

regC2_2 <- feols(Aqua_perc ~ postSurge + postSurge * State + flag | UniqueID + Year, 
                 data = data, vcov = ~UniqueID)
summary(regC2_2)

# apply spatial standard error 
regC2_2c <- feols(Aqua_perc ~ postSurge + postSurge * State + flag | UniqueID + Year, 
                 data = data)
summary(regC2_2c)

se_conley <- vcov_conley(
  regC2_2c ,
  lat = ~Latitude_d,
  lon = ~Longitude_d,
  cutoff = 5  # in kilometers
)
etable(regC2_2c , vcov = se_conley)

# After surge events, aquaculture percentage increases by about 4.217 pp in TN. Significant at 5% 
# In Andhra Pradesh, the surge effect on aquaculture is 3.18 units lower than TN (but not statistically significant).
# In Odisha, the surge effect on aquaculture is 4.26 units lower than TN, which significantly offsets the positive surge effect seen in TN. Significant at 5% 
# => In OD, the surge-driven aquaculture increase seen in TN disappears completely — a statistically significant difference from TN, possibly due to slower recovery, more damage/disruption from surges, Policy or infrastructure differences


# 3. Y = Aquaculture, X = Salinity in previous period (without controls except time and unit, with controls)
regC3_1 <- feols(Aqua_perc ~ avg_salinity_5yr + avg_salinity_5yr * State + flag | Year, 
                 data = data, vcov = ~UniqueID)
summary(regC3_1)

regC3_2 <- feols(Aqua_perc ~ avg_salinity_5yr + avg_salinity_5yr * State + flag | UniqueID + Year, 
                 data = data, vcov = ~UniqueID)
summary(regC3_2)

# apply spatial standard error 
regC3_2c <- feols(Aqua_perc ~ avg_salinity_5yr + avg_salinity_5yr * State + flag | UniqueID + Year, 
                 data = data)
summary(regC3_2c)

se_conley <- vcov_conley(
  regC3_2c ,
  lat = ~Latitude_d,
  lon = ~Longitude_d,
  cutoff = 5  # in kilometers
)
etable(regC3_2c, vcov = se_conley)

# In Tamil Nadu, a 1-unit increase in average salinity over 5 years is associated with a 0.1476 unit decrease in aquaculture percentage, statistically significant at 5% (positive and significant)
# In Andhra Pradesh, the effect of salinity on aquaculture is 0.6493 units higher than in TN (0.6493-0.1476); this difference is significant at 5% (p < 0.05) (positive and significant)
# In Odisha, the effect of salinity on aquaculture is 0.305 units higher than TN, and this difference is significant at 5% (p < 0.05) => +0.1579 (positive and significant)


# 4. Y = Aquaculture, X = Aqua time - 1 (without controls except time and unit, with controls)
regC4_1 <- feols(Aqua_perc ~ Lag_Aqua + Lag_Aqua * State + flag | Year, 
                 data = data, vcov = ~UniqueID)
summary(regC4_1)

regC4_2 <- feols(Aqua_perc ~ Lag_Aqua + Lag_Aqua * State + flag | UniqueID + Year, 
                 data = data, vcov = ~UniqueID)
summary(regC4_2)

# apply spatial standard error 
regC4_2c <- feols(Aqua_perc ~ Lag_Aqua + Lag_Aqua * State + flag | UniqueID + Year, 
                 data = data)
summary(regC4_2c)

se_conley <- vcov_conley(
  regC4_2c ,
  lat = ~Latitude_d,
  lon = ~Longitude_d,
  cutoff = 5  # in kilometers
)
etable(regC4_2c, vcov = se_conley)

# Significant effect of past aquaculture in all states 
# In TN, a 1-unit increase in previous year's aquaculture is associated with a 0.3864 pp increase in current aquaculture. Highly significant.
# In Andhra Pradesh, the effect of past aquaculture is 0.2551 units stronger than in Tamil Nadu. 
# In Odisha, the effect is 0.3014 units stronger than in Tamil Nadu. 



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


models <- list()
models[['Salinity (1)']] <- regC1_2

models[['Aquaculture (1)']] <- regC2_2

models[['Aquaculture (2)']] <- regC3_2

models[['Aquaculture (3)']] <- regC4_2

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
regD1_1 <- feols(Aqua_perc ~ postSurge + avg_salinity_5yr + Lag_Aqua + postSurge * State + flag | Year, 
                 data = data, vcov = ~UniqueID)
summary(regD1_1)

regD1_2 <- feols(Aqua_perc ~ postSurge + avg_salinity_5yr + Lag_Aqua + postSurge * State + flag | UniqueID + Year, 
                 data = data, vcov = ~UniqueID)
summary(regD1_2)

# There is a strong heterogeneity in surge responses across states. 
# Storm surges strongly increase aquaculture in TN, modestly increase it in AP, but no significant surge impact on aquaculture in OD.
# In Tamil Nadu, storm surges are associated with 1.61 point increase in aquaculture share.
# In Andhra Pradesh, Net effect = +0.48
# In Odisha, net effect is = ~0.01 (not statistically different from 0)


# apply spatial standard error 
regD1_2c <- feols(Aqua_perc ~ postSurge + avg_salinity_5yr + Lag_Aqua + postSurge * State + flag | UniqueID + Year, 
                 data = data)
summary(regD1_2c)

se_conley <- vcov_conley(
  regD1_2c ,
  lat = ~Latitude_d,
  lon = ~Longitude_d,
  cutoff = 5  # in kilometers
)
etable(regD1_2c, vcov = se_conley)

# In Tamil Nadu, aquaculture share increases by ~1.6 percentage points after a surge. Significant at 5%.
# No significant direct association between 5-year average salinity and current aquaculture.
# A 1-unit increase in past aquaculture is associated with a 0.6394 unit increase in current aquaculture —> strong persistence.
# In AP, the effect of surge exposure is 1.123 points weaker than in Tamil Nadu, however, this difference is not significant 
# In OD, the effect of surge exposure is 1.591 points weaker than in Tamil Nadu. Although statistically significant at 5%, but the effect is relatively small (0.01)



regD1_3 <- feols(Aqua_perc ~ postSurge + avg_salinity_5yr + Lag_Aqua + avg_salinity_5yr * State + flag | Year, 
                 data = data, vcov = ~UniqueID)
summary(regD1_3)

regD1_4 <- feols(Aqua_perc ~ postSurge + avg_salinity_5yr + Lag_Aqua + avg_salinity_5yr * State + flag | UniqueID + Year, 
                 data = data, vcov = ~UniqueID)
summary(regD1_4)

# Comparing places across the data, in all states higher salinity is related with higher aquaculture, although AP seems to be most responsive to salinity and OD the least. 
# Once unobserved village-level time invarying factors (such as dist or elevation from the sea), etc are controlled for, TN shows a small negative effect of salinity, whereas AP and OD show positive association. 
# In TN, salinity is negatively associated with aquaculture within villages with 0.057 pp
# In AP, salinity is associated with +0.185 pp increase in aquaculture
# In OD, salinity is associated with +0.053 pp increase in aquaculture 


# apply spatial standard error 
regD1_4c <- feols(Aqua_perc ~ postSurge + avg_salinity_5yr + Lag_Aqua + avg_salinity_5yr * State + flag | UniqueID + Year, 
                 data = data)
summary(regD1_4c)

se_conley <- vcov_conley(
  regD1_4c ,
  lat = ~Latitude_d,
  lon = ~Longitude_d,
  cutoff = 5  # in kilometers
)
etable(regD1_4c, vcov = se_conley)

# In Tamil Nadu (reference), storm surge is associated with +0.12pp more aquaculture, but not significant.
# In Tamil Nadu, 5-year average salinity has a negative but association with aquaculture, significant at 5%.
# Strong persistence—past aquaculture drives future aquaculture.
# Salinity’s effect is less negative/more positive in Andhra Pradesh than in TN. 5% significance level. 
# Salinity’s effect is somewhat less negative in Odisha, significant at 10%.



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

print(regD1_2)
print(regD1_4)

# TN aquaculture development seems to be most responsive to shocks, but less so to salinity stress. 
# AP seems to be strategically adaptive to long-term salinity. 
# OD shows least responsiveness to either shocks or stresses, likely due to infrastructure, policy or economic constraints. 

models <- list()
models[['Aquaculture (surge-effect)']] <- regD1_2
models[['Aquaculture (salinity-effect)']] <- regD1_4

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



# Set E : Three-way interactions and non-linearity 
# Try triple interaction: 
regE1_1 <- feols(Aqua_perc ~ postSurge + avg_salinity_5yr + Lag_Aqua + postSurge * avg_salinity_5yr * State + flag | UniqueID + Year, 
                 data = data, vcov = ~UniqueID)
summary(regE1_1)

# Interpretation: 
# TN - Surges have a big positive effect on aquaculture; long-run salinity has a small negative effect, salinity may amplify surge effect, but not significantly. 
# AP - Surges alone are less impacxtful as compared to TN, but salinity is strongly positively associated with aquaculture. Surges in high saline conditions tend to trigger aquaculture expansion the most in AP. 
# OD - Surges have a weak or even negative effect on aqua expansion. Salinity has a more positive response to aquaculture compared to TN. Surges in high saline places have small negative effect, however not significant. 


# Test for non-linear relationship with salinity 
data_clean <- data %>% filter(!is.na(avg_salinity_5yr))
regE1_2 <- feols(Aqua_perc ~ postSurge + poly(avg_salinity_5yr, 2) + Lag_Aqua + postSurge * poly(avg_salinity_5yr, 2) * State + flag | UniqueID + Year, 
                 data = data_clean, vcov = ~UniqueID)
summary(regE1_2)

# Linear effect not significant alone
# Strong curvature effect — suggests increasing returns to salinity for aquaculture in TN.
# Some state-specific nonlinear salinity effects — e.g., salinity effects are stronger OD in higher-order terms.
# A massive, positive coefficient — suggests that in AP, aquaculture after a surge responds very strongly to salinity - maybe visualise it? 

# 5-year average salinity has a nonlinear (convex) effect, meaning aquaculture increases at higher levels of salinity.
# In AP - Aquaculture is highly sensitive to salinity after a surge (strong positive interaction), implying surge and salinity may be joint enablers of expansion.
# In OD - Surge × salinity interactions are not significant. Salinity increases aquaculture moderately, but surge may disrupt or delay response.
# Lagged aquaculture is very predictive — suggesting strong path dependence or persistence in village-level aquaculture.

# apply spatial standard error 
regE1_2c <- feols(Aqua_perc ~ postSurge + poly(avg_salinity_5yr, 2) + Lag_Aqua + postSurge * poly(avg_salinity_5yr, 2) * State + flag | UniqueID + Year, 
                 data = data_clean)
summary(regE1_2c)

se_conley <- vcov_conley(
  regE1_2c ,
  lat = ~Latitude_d,
  lon = ~Longitude_d,
  cutoff = 5  # in kilometers
)
etable(regE1_2c, vcov = se_conley)

# very highly Significant effects: postSurge in TN ; Lag_Aqua ; higher order levels / non-linear salinity 
# Significant at 1%: postsurge impact of salnity in AP
# Significant at 10%: 
# postsurge in TN and OD
# postSurge x poly(...)1 => after surge, salinity has a more positive effect on aquaculture in TN
# StateAP x poly(...)1 => In AP, salinity effect is more positive/nonlinear than in TN
# postSurge x salinity x StateOD => After surge, Odisha sees a reduced or negative nonlinear effect of salinity on aquaculture

# Aquaculture increases after storm surges on average, but the effect is state-specific.
# The relationship between salinity and aquaculture is clearly U-shaped, i.e. Salinity–aquaculture relationship is nonlinear, but varies by state:
# In AP, high salinity after surges strongly increases aquaculture (adaptive behavior).
# In OD, the surge reduces aquaculture (–1.66), with modest nonlinear effects from salinity. The interaction with surge seems to reverse the potential adaptation to aquacultre - presents itself as a shock. 
# No strong national-level salinity × surge interaction, but clear heterogeneity across states.

# Policy implications
# Supports the hypothesis of adaptive aquaculture expansion in response to salinization, especially in AP after surges.
# Emphasizes the importance of regional differences — Odisha may not have similar coping strategies or infrastructure.
# Suggests salinity impacts are nonlinear, which is important for designing thresholds/interventions.


models <- list()
models[['Aquaculture (three-way interaction)']] <- regE1_1
models[['Aquaculture (non-linearity)']] <- regE1_2

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


# If need to restrict to villages only with aqua 
villages_with_aqua <- data %>%
  group_by(UniqueID) %>%
  summarize(has_aqua = any(!is.na(Aqua_perc) & Aqua_perc > 0)) %>%
  filter(has_aqua) %>%
  pull(UniqueID)

data_filtered <- data %>%
  filter(UniqueID %in% villages_with_aqua)



# 1. Y = Current Salinity, X = Lag Aqua (without controls except time and unit, with controls)
regF1_1 <- feols(Saline_perc_norm  ~ Lag_Aqua | UniqueID + Year, 
                 data = data, vcov = ~UniqueID)
summary(regF1_1)

regF1_2 <- feols(Saline_perc_norm  ~ Lag_Aqua + postSurge | UniqueID + Year, 
                 data = data, vcov = ~UniqueID)
summary(regF1_2)
# Adding postSurge has no impact on the model or the coefficient of past aquaculture 

# apply spatial standard error 
regF1_2c <- feols(Saline_perc_norm  ~ Lag_Aqua + postSurge | UniqueID + Year, 
                 data = data)
summary(regF1_2c)

se_conley <- vcov_conley(
  regF1_2c,
  lat = ~Latitude_d,
  lon = ~Longitude_d,
  cutoff = 5  # in kilometers
)
etable(regF1_2c, vcov = se_conley)

# Effects still significant, although at a lower significance levels

models <- list()
models[['Salinity (1)']] <- regF1_1
models[['Salinity (2)']] <- regF1_2
msummary(models, stars = c('*' = .1, '**' = .05, '***' = .01),gof_omit=c("BIC|AIC|RMSE|R2 Within Adj."),coef_omit=c("(Intercept)"), filename = 'table.rtf')




# 2. Y = Future avg salinity, X = Current Aqua
regF2_1 <- feols(future_avg_salinity ~ Aqua_perc | UniqueID + Year, 
                 vcov = ~UniqueID, data = data)

summary(regF2_1)

regF2_2 <- feols(future_avg_salinity ~ Aqua_perc + postSurge | UniqueID + Year, 
                 vcov = ~UniqueID, data = data)

summary(regF2_2)

models <- list()
models[['Future Salinity (1)']] <- regF2_1
models[['Future Salinity (2)']] <- regF2_2
msummary(models, stars = c('*' = .1, '**' = .05, '***' = .01),gof_omit=c("BIC|AIC|RMSE|R2 Within Adj."),coef_omit=c("(Intercept)"), filename = 'table.rtf')


# apply spatial standard error 
regF2_2c <- feols(future_avg_salinity ~ Aqua_perc + postSurge | UniqueID + Year, 
                data = data)

summary(regF2_2c)

se_conley <- vcov_conley(
  regF2_2c,
  lat = ~Latitude_d,
  lon = ~Longitude_d,
  cutoff = 5  # in kilometers
)
etable(regF2_2c, vcov = se_conley)
# Past aquaculture is a significant predictor of future salinity, despite severe controls and spatial error handling



# These results provide strong empirical evidence of a lagged ecological feedback from human adaptation (aquaculture expansion) to system degradation (salinization). Even after removing the "zero-aquaculture" villages (which might dilute the effect), you still observe a statistically robust and consistent association.
# The effect persists across model specifications; The postSurge control confirms the model isn't confounded by extreme weather shocks
# The magnitude is modest but credible for environmental change processes. It is likley modest, because the impact of aquculture on salinity is more localised (in the absolute neighborhoods of the aquaculture ponds which is averaged across the village)


# 3. By State
regF3_1 <- feols(Saline_perc_norm  ~ Lag_Aqua * State | UniqueID + Year, 
                 data = data, vcov = ~UniqueID)
summary(regF3_1)

regF3_2 <- feols(Saline_perc_norm  ~ Lag_Aqua * State + postSurge | UniqueID + Year, 
                 data = data, vcov = ~UniqueID)
summary(regF3_2)

# Adding postSurge has no impact on the model explanation or the coefficient of past aquaculture 

# apply spatial standard error 
regF3_2c <- feols(Saline_perc_norm  ~ Lag_Aqua * State + postSurge | UniqueID + Year, 
                 data = data)
summary(regF3_2c)

se_conley <- vcov_conley(
  regF3_2c,
  lat = ~Latitude_d,
  lon = ~Longitude_d,
  cutoff = 5  # in kilometers
)
etable(regF3_2c, vcov = se_conley)

# Higher past aquaculture slightly reduces salinity in the next period in TN. Significant at 1%
# After a surge, normalized salinity drops slightly (could indicate flushing, dilution, or temporary land abandonment). Significant at 5% 
# In AP and OD, past aquaculture is positively associated with future salinity – indicating possible degradation or intensification effects. Highly significant in both


models <- list()
models[['Salinity (1)']] <- regF3_1
models[['Salinity (2)']] <- regF3_2
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
a$Variable <- factor(a$Variable, levels=c("Aquaculture (Past)","Surge and Salinity","Salinity (in non-surge areas)","Surge (in low saline areas)"))
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
