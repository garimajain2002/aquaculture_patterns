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


# Apart from Tiem and Place FE, we also need to control for place-specific time variant factors that may have their own impact on aquaculture.
# These could be state level policies in a year, level of urbanization and infrastructure, or rainfall (more rains => better agricultural output = less likely to shift) 
# Add these as additional controls and check results


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


# Convert Population to Density first 
data$density <- (data$Population / data$Shape_Area) * 10000
summary(data$density)

# Assign district level rainfall averages to the few misisng rainfall points 
summary(data$Rainfall_mm)
data <- data %>%
  group_by(District) %>%
  mutate(
    district_mean_rainfall = mean(Rainfall_mm, na.rm = TRUE),
    Rainfall_mm = if_else(is.na(Rainfall_mm), district_mean_rainfall, Rainfall_mm)
  ) %>%
  ungroup() %>%
  select(-district_mean_rainfall)  # remove helper column

summary(data$Rainfall_mm)

# Convert rainfall from mm to m to make the results more meaningful
data$Rainfall_m <- data$Rainfall_mm / 1000


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

# with controls
# Rainfall and density might also affect the level of salinity. Can be added as control 
regA1_2c <- feols(avg_salinity_5yr ~ postSurge + flag + Rainfall_m + density | UniqueID + Year, 
                 data = data, vcov = ~UniqueID)
summary(regA1_2c)

# with state specific time variation 
regA1_2cS <- feols(avg_salinity_5yr ~ postSurge + flag + Rainfall_m + density + i(State, Year, ref = "TN") | UniqueID + Year, 
                  data = data, vcov = ~UniqueID)
summary(regA1_2cS)


# The comparison between without and with state-level variation highlights that post-surge variable tends to capture state variations in post-surge responses that might affect aquaculture uptake. 
# post-surge is therefore a variable that depicts not just the surge effect, but also a shift in how the states respond and thereby affecting the development trajectories of these places. 
# Remember, not each surge is the same in terms of its impact. Diufferent levels of disasters attract different actions from state and non-state actors 


# apply spatial standard error 

regA1_2cc <- feols(avg_salinity_5yr ~ postSurge + flag + Rainfall_m + density | UniqueID + Year, 
                 data = data)
summary(regA1_2cc)

se_conley <- vcov_conley(
  regA1_2cc,
  lat = ~Latitude_d,
  lon = ~Longitude_d,
  cutoff = 5  # in kilometers
)
etable(regA1_2cc, vcov = se_conley)

# The coefficient of postSurge is not statistically significant, ie, the ratio of the coefficient to its SE (t-statistic) is less than expected (1.65 for 10%, 1.96 for 5%, and 2.58 for 1% significance) 
# This is mainly because the effect is too small wrt to the SE 
# Therefore, we cannot reject the null hypothesis that postSurge has no effect on avg_salinity_5yr.
# See geographical heterogeneity in this - AP is significant at 1%, TN at 5% and OD not significantly different from TN. 
# Density is the most statistically significant predictor of salinity (i.e. higher urbanization, higher salinity) 

# apply spatial error with state year variations
regA1_2ccS <- feols(avg_salinity_5yr ~ postSurge + flag + Rainfall_m + density + i(State, Year, ref ="TN") | UniqueID + Year, 
                    data = data)
summary(regA1_2ccS)

se_conley <- vcov_conley(
  regA1_2ccS,
  lat = ~Latitude_d,
  lon = ~Longitude_d,
  cutoff = 5  # in kilometers
)
etable(regA1_2ccS, vcov = se_conley)

# Odisha has on average 0.0064 units lower salinity than Tamil Nadu over time (net of year and place FEs). This is statistically significant, suggesting structural differences in baseline salinity across states.



# Similarly, effect of rainfall on salinity 
# 1R. Y = Salinity, X = surge (without controls except time and unit)
regA1_1_1 <- feols(avg_salinity_5yr ~ Rainfall_m + flag | Year, 
                 data = data, vcov = ~UniqueID)
summary(regA1_1_1)

regA1_2_1 <- feols(avg_salinity_5yr ~ Rainfall_m + flag | UniqueID + Year, 
                 data = data, vcov = ~UniqueID)
summary(regA1_2_1)

# Time FE: 
# Rainfall seems to be significant predictor 

# Village + Time FE
# Significance goes away => 

# The flip suggests that areas that get surges more often may already have higher salinity. The negative effect may indicate a temporary leaching/flushing effect or behavioral change (e.g., fallowing land post-surge) or higher soil moisture

# with state controls

regA1_2_1cS <- feols(avg_salinity_5yr ~ Rainfall_m + flag + i(State, Year, ref = "TN") | UniqueID + Year, 
                   data = data, vcov = ~UniqueID)
summary(regA1_2_1cS)



# apply spatial standard error 

regA1_2_1ccS <- feols(avg_salinity_5yr ~ Rainfall_m + flag + i(State, Year, ref = "TN") | UniqueID + Year, 
                    data = data)
summary(regA1_2_1ccS)


se_conley <- vcov_conley(
  regA1_2_1ccS,
  lat = ~Latitude_d,
  lon = ~Longitude_d,
  cutoff = 5  # in kilometers
)
etable(regA1_2_1ccS, vcov = se_conley)

# The coefficient of rainfall is not statistically significant => has no association with levels of salinity 
# However the baseline salinity in Odisha is lower from TN, and thsi difference is statistically different




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


# with controls 
regA2_2c <- feols(Aqua_perc ~ postSurge + flag + Rainfall_m + density | UniqueID + Year, 
                 data = data, vcov = ~UniqueID)
summary(regA2_2c)

# with State controls 
regA2_2cS <- feols(Aqua_perc ~ postSurge + flag + Rainfall_m + density + i(State, Year, ref="TN") | UniqueID + Year, 
                  data = data, vcov = ~UniqueID)
summary(regA2_2cS)

# It does not change the model fit very much, although takes away the significance of postSurge => postSurge in teh previous model stands for state variations in responses



# apply spatial standard error 
regA2_2cc <- feols(Aqua_perc ~ postSurge + flag + Rainfall_m + density | UniqueID + Year, 
                 data = data)
summary(regA2_2cc)

se_conley <- vcov_conley(
  regA2_2cc,
  lat = ~Latitude_d,
  lon = ~Longitude_d,
  cutoff = 5  # in kilometers
)
etable(regA2_2cc, vcov = se_conley)

# The ratio with SE (0.1889) is ~1.221 which is below the 1.645 for 10% significance 
# Cannot exclude the null hypothesis that surge has no effect on aquaculture uptake once spatial clustering is accounted for


regA2_2ccS <- feols(Aqua_perc ~ postSurge + flag + Rainfall_m + density + i(State, Year, ref="TN")| UniqueID + Year, 
                   data = data)
summary(regA2_2ccS)

se_conley <- vcov_conley(
  regA2_2ccS,
  lat = ~Latitude_d,
  lon = ~Longitude_d,
  cutoff = 5  # in kilometers
)
etable(regA2_2ccS, vcov = se_conley)

# density is a predictor of how aquaculture is adopted 
# State variations in aquaculture adoption baseline - For Odisha, aquaculture share is .04 pp higher, on average, in each year than TN.




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


# with controls 
regA3_2c <- feols(Aqua_perc ~ avg_salinity_5yr + flag + Rainfall_m + density | UniqueID + Year, 
                 data = data, vcov = ~UniqueID)
summary(regA3_2c)

# with state time variations controlled 
regA3_2cS <- feols(Aqua_perc ~ avg_salinity_5yr + flag + Rainfall_m + density + i(State, Year, ref = "TN") | UniqueID + Year, 
                  data = data, vcov = ~UniqueID)
summary(regA3_2cS)


# apply spatial standard error 
regA3_2cc <- feols(Aqua_perc ~ avg_salinity_5yr + flag + Rainfall_m + density | UniqueID + Year, 
                 data = data)
summary(regA3_2cc)

se_conley <- vcov_conley(
  regA3_2cc,
  lat = ~Latitude_d,
  lon = ~Longitude_d,
  cutoff = 5  # in kilometers
)
etable(regA3_2cc, vcov = se_conley)

# A 10 unit increase in salinity over the past 5 years is associated with a 2 pp increase in aquaculture area. 
# This is significant (p < 0.1), perhaps because of a potential long-run feedback loop from salinity to aquaculture expansion.
# Rainfall does not seem to play a meaningful role here — this might be because the effect is either truly null, or already captured by other spatial/temporal variation.
# Higher population density is associated with less aquaculture, which makes sense in more urbanized/dense areas.

regA3_2ccS <- feols(Aqua_perc ~ avg_salinity_5yr + flag + Rainfall_m + density + i(State, Year, ref="TN") | UniqueID + Year, 
                   data = data)
summary(regA3_2ccS)

se_conley <- vcov_conley(
  regA3_2ccS,
  lat = ~Latitude_d,
  lon = ~Longitude_d,
  cutoff = 5  # in kilometers
)
etable(regA3_2ccS, vcov = se_conley)




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


# with controls
regA4_2c <- feols(Aqua_perc ~ Lag_Aqua + flag + Rainfall_m + density | UniqueID + Year, 
                 data = data, vcov = ~UniqueID)
summary(regA4_2c)


# with state time variation controlled 
regA4_2cS <- feols(Aqua_perc ~ Lag_Aqua + flag + Rainfall_m + density + i(State, Year, ref = "TN") | UniqueID + Year, 
                  data = data, vcov = ~UniqueID)
summary(regA4_2cS)



# apply spatial standard error 
regA4_2cc <- feols(Aqua_perc ~ Lag_Aqua + flag + Rainfall_m + density | UniqueID + Year, 
                 data = data)
summary(regA4_2cc)

se_conley <- vcov_conley(
  regA4_2cc,
  lat = ~Latitude_d,
  lon = ~Longitude_d,
  cutoff = 5  # in kilometers
)
etable(regA4_2cc, vcov = se_conley)


# The estimate for lag_aquua is highly statistically significant. 
# strong evidence of path dependence or inertia in aquaculture. Once adopted, aquaculture practices tend to persist, possibly due to:
# Sunk investments in ponds or infrastructure, Skill specialization or market commitments, Environmental lock-ins (e.g., salinized land).
# Rainfall is not not sigificant 
# Density remains negatively associated, but the magnitude is small

# apply spatial standard error 
regA4_2ccS <- feols(Aqua_perc ~ Lag_Aqua + flag + Rainfall_m + density + i(State, Year, ref = "TN") | UniqueID + Year, 
                   data = data)
summary(regA4_2ccS)

se_conley <- vcov_conley(
  regA4_2ccS,
  lat = ~Latitude_d,
  lon = ~Longitude_d,
  cutoff = 5  # in kilometers
)
etable(regA4_2ccS, vcov = se_conley)




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

# with controls 
regA5_2c <- feols(Aqua_perc ~ postSurge + avg_salinity_5yr + Lag_Aqua + flag + Rainfall_m + density | UniqueID + Year, 
                 data = data, vcov = ~UniqueID)
summary(regA5_2c)


# with state time variation 
regA5_2cS <- feols(Aqua_perc ~ postSurge + avg_salinity_5yr + Lag_Aqua + flag + Rainfall_m + density + i(State, Year, ref = "TN") | UniqueID + Year, 
                  data = data, vcov = ~UniqueID)
summary(regA5_2cS)



# apply spatial standard error 
regA5_2cc <- feols(Aqua_perc ~ postSurge + avg_salinity_5yr + Lag_Aqua + flag + Rainfall_m + density | UniqueID + Year, 
                 data = data)
summary(regA5_2cc)

se_conley <- vcov_conley(
  regA5_2cc,
  lat = ~Latitude_d,
  lon = ~Longitude_d,
  cutoff = 5  # in kilometers
)
etable(regA5_2cc, vcov = se_conley)

# Slight increase in aquaculture after surge, weakly significant.
# Higher salinity over past 5 years is weakly associated with more aquaculture.
# Very strong and highly significant autocorrelation — villages with aquaculture in t−1 tend to continue it in t.

regA5_2ccS <- feols(Aqua_perc ~ postSurge + avg_salinity_5yr + Lag_Aqua + flag + Rainfall_m + density + i(State, Year, ref = "TN") | UniqueID + Year, 
                   data = data)
summary(regA5_2ccS)

se_conley <- vcov_conley(
  regA5_2ccS,
  lat = ~Latitude_d,
  lon = ~Longitude_d,
  cutoff = 5  # in kilometers
)
etable(regA5_2ccS, vcov = se_conley)


 
models <- list()
models[['Salinity 1 (Time FE)']] <- regA1_1
models[['Salinity 1 (Village+Time FE)']] <- regA1_2
models[['Salinity 1 (Village+Time FE+Controls)']] <- regA1_2c

models[['AC 1 (Time FE)']] <- regA2_1
models[['AC 1 (Village+Time FE)']] <- regA2_2
models[['AC 1 (Village+Time FE + Controls)']] <- regA2_2c

models[['AC 2 (Time FE)']] <- regA3_1
models[['AC 2 (Village+Time FE)']] <- regA3_2
models[['AC 2 (Village+Time FE + Controls)']] <- regA3_2c

models[['AC 3 (Time FE)']] <- regA4_1
models[['AC 3 (Village+Time FE)']] <- regA4_2
models[['AC 3 (Village+Time FE + Controls)']] <- regA4_2c

models[['AC 4 (Time FE)']] <- regA5_1
models[['AC 4 (Village+Time FE)']] <- regA5_2
models[['AC 4 (Village+Time FE + Controls)']] <- regA5_2c

msummary(models, stars = c('*' = .1, '**' = .05, '***' = .01),gof_omit=c("BIC|AIC|RMSE|R2 Within Adj."),coef_omit=c("(Intercept)"), filename = 'table.rtf')


# Adding the controls doesn't improve the model fit or the other coefficients very much 



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

# with controls
regB1_2c <- feols(Aqua_perc ~ postSurge + avg_salinity_5yr + Lag_Aqua + postSurge * avg_salinity_5yr + flag  + Rainfall_m + density | UniqueID + Year, 
                 data = data, vcov = ~UniqueID)
summary(regB1_2c)

# with state time controls
regB1_2cS <- feols(Aqua_perc ~ postSurge + avg_salinity_5yr + Lag_Aqua + postSurge * avg_salinity_5yr + flag  + Rainfall_m + density + i(State, Year, ref = "TN") | UniqueID + Year, 
                  data = data, vcov = ~UniqueID)
summary(regB1_2cS)




# apply spatial standard error 
regB1_2cc <- feols(Aqua_perc ~ postSurge + avg_salinity_5yr + Lag_Aqua + postSurge * avg_salinity_5yr + flag + Rainfall_m + density | UniqueID + Year, 
                 data = data)
summary(regB1_2cc)

se_conley <- vcov_conley(
  regB1_2cc,
  lat = ~Latitude_d,
  lon = ~Longitude_d,
  cutoff = 5  # in kilometers
)
etable(regB1_2cc, vcov = se_conley)

# Only lag aqua and density are significant. But also, within R2 suggests that 38% of within village variation is captured by these measures. Adding rainfall and density doesnt change the model fit. 


regB1_2ccS <- feols(Aqua_perc ~ postSurge + avg_salinity_5yr + Lag_Aqua + postSurge * avg_salinity_5yr + flag + Rainfall_m + density + i(State, Year, ref = "TN") | UniqueID + Year, 
                   data = data)
summary(regB1_2ccS)

se_conley <- vcov_conley(
  regB1_2ccS,
  lat = ~Latitude_d,
  lon = ~Longitude_d,
  cutoff = 5  # in kilometers
)
etable(regB1_2ccS, vcov = se_conley)

# Only lag aqua and density are significant. But also, within R2 suggests that 38% of within village variation is captured by these measures. Adding rainfall and density doesnt change the model fit. 


models <- list()
models[['AC (All controls+Time FE)']] <- regB1_1
models[['AC (All controls+Village+Time FE)']] <- regB1_2

msummary(models, stars = c('*' = .1, '**' = .05, '***' = .01),gof_omit=c("BIC|AIC|RMSE|R2 Within Adj."),coef_omit=c("(Intercept)"), filename = 'table.rtf')


# without place-time varying controls 
models <- list()
models[['Salinity (1)']] <- regA1_2

models[['Aquaculture (1)']] <- regA3_2

models[['Aquaculture (2)']] <- regA2_2

models[['Aquaculture (3)']] <- regA4_2

models[['Aquaculture (4-All)']] <- regA5_2

models[['Aquaculture (5-Interaction)']] <- regB1_2

msummary(models, stars = c('*' = .1, '**' = .05, '***' = .01),gof_omit=c("BIC|AIC|RMSE|R2 Within Adj."),coef_omit=c("(Intercept)"), filename = 'table.rtf')



# with additional place-time varying controls 
models <- list()
models[['Salinity (1)']] <- regA1_2c

models[['Aquaculture (1)']] <- regA3_2c

models[['Aquaculture (2)']] <- regA2_2c

models[['Aquaculture (3)']] <- regA4_2c

models[['Aquaculture (4-All)']] <- regA5_2c

models[['Aquaculture (5-Interaction)']] <- regB1_2c

msummary(models, 
         stars = c('*' = .1, '**' = .05, '***' = .01),
         gof_omit = c("BIC|AIC|RMSE|R2 Within Adj."),
         coef_omit = c("(Intercept)"), 
         coef_map = c(
           "postSurge" = "Post-Surge",
           "avg_salinity_5yr" = "Persistent Salinity",
           "postSurge:avg_salinity_5yr" = "Surge in High Saline Areas",
           "Rainfall_m" = "Rainfall",
           "density" = "Population Density",
           "Lag_Aqua" = "Aquaculture (Previous Year)"
         ),
         filename = 'table.rtf'

)


# with state time varying controls 
models <- list()
models[['Salinity (1)']] <- regA1_2cS

models[['Aquaculture (1)']] <- regA3_2cS

models[['Aquaculture (2)']] <- regA2_2cS

models[['Aquaculture (3)']] <- regA4_2cS

models[['Aquaculture (4-All)']] <- regA5_2cS

models[['Aquaculture (5-Interaction)']] <- regB1_2cS

msummary(models, 
         stars = c('*' = .1, '**' = .05, '***' = .01),
         gof_omit = c("BIC|AIC|RMSE|R2 Within Adj."),
         coef_omit = c("(Intercept)"), 
         coef_map = c(
           "postSurge" = "Post-Surge",
           "avg_salinity_5yr" = "Persistent Salinity",
           "postSurge:avg_salinity_5yr" = "Surge in High Saline Areas",
           "Rainfall_m" = "Rainfall",
           "density" = "Population Density",
           "Lag_Aqua" = "Aquaculture (Previous Year)"
         ),
         filename = 'table.rtf'
         
)




# Print selected models for comparison - with increasing complexity and spatial robustness 
# Predictors + FE with Clustered SE
# Predictors + Interaction + FE with Clustered SE
# Predictors + Interaction + Controls + FE with Clustered SE
# Predictors + Interaction + Controls + FE with Conley SE 

regA5_2 <- feols(Aqua_perc ~ postSurge + avg_salinity_5yr + Lag_Aqua + flag | UniqueID + Year, 
                 data = data, vcov = ~UniqueID)
summary(regA5_2)

regB1_2 <- feols(Aqua_perc ~ postSurge + avg_salinity_5yr + Lag_Aqua + postSurge * avg_salinity_5yr + flag | UniqueID + Year, 
                 data = data, vcov = ~UniqueID)
summary(regB1_2)

regB1_2c <- feols(Aqua_perc ~ postSurge + avg_salinity_5yr + Lag_Aqua + postSurge * avg_salinity_5yr + flag  + Rainfall_m + density | UniqueID + Year, 
                  data = data, vcov = ~UniqueID)
summary(regB1_2c)

regB1_2cS <- feols(Aqua_perc ~ postSurge + avg_salinity_5yr + Lag_Aqua + postSurge * avg_salinity_5yr + flag  + Rainfall_m + density + i(State, Year, ref = "TN") | UniqueID + Year, 
                   data = data, vcov = ~UniqueID)
summary(regB1_2cS)

regB1_2ccS <- feols(Aqua_perc ~ postSurge + avg_salinity_5yr + Lag_Aqua + postSurge * avg_salinity_5yr + flag + Rainfall_m + density + i(State, Year, ref = "TN") | UniqueID + Year, 
                    data = data)
summary(regB1_2ccS)

se_conley <- vcov_conley(
  regB1_2ccS,
  lat = ~Latitude_d,
  lon = ~Longitude_d,
  cutoff = 5  # in kilometers
)
etable(regB1_2ccS, vcov = se_conley)


models <- list(
  "Aquaculture (1) Predictors + TWFE" = regA5_2,
  "Aquaculture (2) Predictors + Interaction + TWFE" = regB1_2,
  "Aquaculture (3) Predictors + Interaction + Controls + TWFE" = regB1_2c,
  "Aquaculture (4) Predictors + Interaction + Controls + State-Time factors + TWFE" = regB1_2cS,
  "Aquaculture (5) Predictors + Interaction + Controls + State-Time factors + TWFE (Spatial)" = regB1_2ccS
)

se_list <- list(
  "Clustered" = ~UniqueID, 
  "Clustered" = ~UniqueID,
  "Clustered" = ~UniqueID,
  "Clustered" = ~UniqueID,
  "Spatial" = vcov_conley(
    regB1_2ccS,
    lat = ~Latitude_d,
    lon = ~Longitude_d,
    cutoff = 5
  )
)

msummary(
  models,
  vcov = se_list,
  gof_omit = "Adj|Within|Log|AIC|BIC", 
  stars = c('*' = 0.1, '**' = 0.05, '***' = 0.01)
)


# Comparing Model 3 and 4 -
# Adding State × Year fixed effects in Model 4 reveals that the previously positive average association between storm surges and aquaculture disappears — even turning negative — once state-specific annual shocks are controlled for. 
# However, the interaction between storm surges and persistent salinity becomes even stronger, suggesting that aquaculture expansion is not a generalized response to storm events, but a targeted adaptation in areas where salinity stress is already high. 
# The stable model fit (R²) across specifications indicates that this interaction captures the key localized environmental mechanism driving change, beyond what is explained by broader policy or macroeconomic trends.


# The fact that Conley SEs are larger means that spatial correlation is present — nearby villages are likely experiencing similar shocks or environmental conditions.
# In other words, storm surges and salinity cluster spatially. This validates the core hypothesis that aquaculture expansion is geographically concentrated in areas affected by surges and salinity.
# Significance disappears not because the effect is gone, but because standard errors are now more realistic — accounting for spatial dependence that clustered SEs alone miss.
# The effect sizes (coefficients) remain nearly identical even after controlling for State × Year and spatial SEs — they don't shrink.
# The R² stays flat (~0.746–0.747) across all models. This shows: the model fit is stable. Adding spatial controls doesn’t explain much more variation, meaning that postSurge and salinity already capture meaningful variation in aquaculture.

# Although the statistical significance of postSurge and persistent salinity weakens under spatial standard errors, 
# the coefficient magnitudes remain robust across specifications. This suggests that their impact on aquaculture is not spurious but embedded in spatially clustered dynamics. 
# In fact, the loss of precision reflects the very nature of the processes we study: storm surges and salinization are spatially concentrated phenomena. 
# The consistency in R² and effect sizes across models confirms the substantive importance of these variables.






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


models <- list()
models[['Salinity (1)']] <- regF1_2
models[['Salinity (2)']] <- regF3_2
msummary(models, stars = c('*' = .1, '**' = .05, '***' = .01),gof_omit=c("BIC|AIC|RMSE|R2 Within Adj."),coef_omit=c("(Intercept)"), filename = 'table.rtf')




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


