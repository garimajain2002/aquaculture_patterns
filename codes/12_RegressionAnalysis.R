# Aquaculture Analysis: Panel Data Regression Models
# This script analyzes the relationship between aquaculture land change, storm surges, salinity, and other socio-environmental factors

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

getwd()

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

# Create 5-year moving average for salinity if needed
data <- data %>%
  mutate(Year = as.numeric(Year)) %>%  # Ensure Year is numeric
  arrange(UniqueID, Year) %>%
  group_by(UniqueID) %>%
  complete(Year = full_seq(Year, 1)) %>%  # Ensure continuous years for each village
  mutate(
    avg_salinity_5yr = rollapply(Saline_perc_norm, 
                                 width = 5, 
                                 FUN = mean, 
                                 fill = NA, 
                                 align = "right", 
                                 na.rm = TRUE, 
                                 partial = TRUE)  # Allow smaller windows when data is missing
  ) %>%
  ungroup()

summary(data$avg_salinity_5yr)
summary(data$Avg_Salinity_Last5)
# Check how is this diff from the previously calculated Avg_Salinity_Last5
# fewer NAs in this version. Rest the same. Use "avg_salinity_5yr" for analysis 


# Create interaction term for ease
data$salinity_storm_interaction <- data$Saline_perc_norm * data$postSurge
summary(data)


# ----------------------
# REGRESSION MODELS
# ----------------------
# RQ: In what ways is the area under aquaculture related to salinity prevalence and history of storm surges,
# accounting for distance from the sea and across-state variations. 

# Model 1: Basic pooled OLS regression
model_1 <- lm(Aqua_perc ~ postSurge + Saline_perc_norm + 
               Sea_Dist + DEM_avg, data = data)
summary(model_1)

# Model 2: Adding the interaction term
model_2 <- lm(Aqua_perc ~ postSurge + Saline_perc_norm + 
               salinity_storm_interaction + Sea_Dist + DEM_avg, data = data)
summary(model_2)


# Model 3 - State level FE + Two way clustering (also clustering for villages across years)
# Relevel State so that TN is the reference category
data$State <- relevel(factor(data$State), ref = "TN")

model_3 <- feols(Aqua_perc ~ postSurge + Saline_perc_norm + 
                   salinity_storm_interaction + (State * postSurge) + DEM_avg | Sea_Dist, 
                 data = data, 
                 vcov = ~UniqueID ) 
summary(model_3)


models <- list()
models[['Basic OLS']] <- model_1
models[['OLS + Interaction']] <- model_2
models[['Fixed Effects (Dist) + Clustered SE']] <- model_3
msummary(models, stars = c('*' = .1, '**' = .05, '***' = .01),gof_omit=c("BIC|AIC|RMSE|R2 Within Adj."),coef_omit=c("(Intercept)"), filename = 'table.rtf')

# Interpretation: 
# In TN, storm surges are associated with a +0.53 pp increase in aquaculture %, at average salinity and elevation.
# A 1 SD increase in normalized salinity is associated with a +0.91 pp increase in aquaculture %.
# The effect of salinity on aquaculture is less positive after a storm surge—by ~0.463 pp per SD.
# In non-surge years, Andhra Pradesh has +0.66 pp higher aquaculture % than TN.
# In non-surge years, Odisha has +0.21 pp higher aquaculture % than TN.
# The storm surge effect in AP is 4.04 pp greater than in TN. Total effect = 0.53 + 4.04 ≈ +4.59 pp
# The storm surge effect in OD is 1.03 pp less than in TN. Total effect = 0.524 − 1.03 ≈ −0.51 pp
# Aquaculture slightly declines with elevation — by 0.06 pp for every 100 meters of elevation.
# RMSE: ~5.2, meaning the average model prediction error is ~5 percentage points.
# Adjusted R²: ~0.051 — low, but reasonable given village-level variation.
# Within R²: ~0.048 — ~4.8% of within-village variation in aquaculture is explained.

# This suggests that storm surges are strongly associated with increased aquaculture in Andhra Pradesh, modest increases in Tamil Nadu, and even declines in Odisha. 
# Higher salinity levels also generally promote aquaculture, but their positive impact weakens after a surge. 
# Slight negative effects of elevation suggest topographic barriers to expansion. 
# Overall, there are visible state-level differences in how aquaculture responds to climate shocks, shaped by both geography, policies, and prior land use dynamics.


# For persistent Salinity - Using the 5-year average salinity
# Model 4: Basic pooled OLS regression
model_4 <- lm(Aqua_perc ~ postSurge + avg_salinity_5yr + 
                Sea_Dist + DEM_avg, data = data)
summary(model_4)

# Model 5: Adding the interaction term
model_5 <- lm(Aqua_perc ~ postSurge + avg_salinity_5yr + 
                avg_salinity_5yr*postSurge + Sea_Dist + 
                DEM_avg, data = data)
summary(model_5)


# Model 6: State level variations + clustering for villages across years
model_6 <- feols(Aqua_perc ~ postSurge + avg_salinity_5yr + 
                   (avg_salinity_5yr*postSurge) + (State * postSurge) + DEM_avg | Sea_Dist, 
                 data = data, 
                 vcov = ~UniqueID )
summary(model_6)

models <- list()
models[['Basic OLS']] <- model_4
models[['OLS + Interaction']] <- model_5
models[['Fixed Effects (Dist) + Clustered SE']] <- model_6
msummary(models, stars = c('*' = .1, '**' = .05, '***' = .01),gof_omit=c("BIC|AIC|RMSE|R2 Within Adj."),coef_omit=c("(Intercept)"), filename = 'table.rtf')
# Interpretation: 
# Before a storm surge, a 1 unit increase in average salinity over the last 5 years is associated with a 1.46 percentage point increase in aquaculture.
# suggests strong association of prior salinity levels with aquaculture.
# After a storm surge, aquaculture increased by 0.47 pp on average in TN, holding average salinity constant.
# StateAP: +0.58* → Andhra Pradesh villages have higher baseline aquaculture than TN, even before the storm.
# StateOD: +0.12 **→ Odisha also shows higher aquaculture than TN pre-storm, but much smaller effect.
# postSurge:StateAP: +4.03* → AP sees a much larger increase (4.03 pp) in aquaculture post-surge compared to TN.
# postSurge:StateOD: –1.15* → Odisha sees a decline of ~1.15 pp relative to TN in post-surge aquaculture.
# Interaction of salinity and storm : The effect of average salinity on aquaculture declines by 0.78 pp after a storm.
# So post-surge, the marginal effect of salinity becomes: 1.46 – 0.78 ≈ 0.68 pp increase per unit salinity → salinity still increases aquaculture, but less strongly.
# Aquaculture is strongly predicted by long-term salinity, but that relationship weakens after storm surges. 
# Andhra Pradesh shows much greater post-surge aquaculture expansion than Tamil Nadu, while Odisha shows a drop. 
# Salinity, even averaged over 5 years, encourages aquaculture, but its effect is tempered after surge events, possibly due to damage, altered soil dynamics, or adaptation lag.



# With Salinity as binary 
# Model 7: Basic pooled OLS regression
model_7 <- lm(Aqua_perc ~ postSurge + Saline + 
                Sea_Dist + DEM_avg, data = data)
summary(model_7)

# Model 8: Adding the interaction term
model_8 <- lm(Aqua_perc ~ postSurge + Saline + 
                Saline * postSurge + Sea_Dist+ 
                DEM_avg, data = data)
summary(model_8)


# Model 9 - State level variations + FE for distance + clustering errors for villages across years
model_9 <- feols(Aqua_perc ~ postSurge + Saline + 
                   (Saline*postSurge) + (State * postSurge) + DEM_avg | Sea_Dist, 
                 data = data, 
                 vcov = ~UniqueID )
summary(model_9)

models <- list()
models[['Basic OLS']] <- model_7
models[['OLS + Interaction']] <- model_8
models[['Fixed Effects (State+Dist)']] <- model_9

msummary(models, stars = c('*' = .1, '**' = .05, '***' = .01),gof_omit=c("BIC|AIC|RMSE|R2 Within Adj."),coef_omit=c("(Intercept)"), filename = 'table.rtf')
# Interpretation: In pre-surge periods, villages above the median salinity level have 0.57 percentage points more aquaculture than low-salinity villages.
# PostSurge coefficient is for low-salinity village, i.e. In low-salinity villages, the storm surge is associated with a +0.20 pp increase in aquaculture. 
# This is not statistically significant, implies not much evidence of change in low-saline villages after a surge.
# Interaction belween postSurge and Saline shows the additional effect of a storm surge in high-saline villages.
# High-saline villages see an extra 0.55 pp increase in aquaculture after a storm compared to low-saline villages.
# post-surge aquaculture in high-saline villages rises by: postSurge + postSurge:Saline = 0.20 + 0.55 = ~0.75 pp
# Andhra Pradesh and Odisha both have higher baseline aquaculture (pre-surge) than Tamil Nadu (reference state), especially AP.
# Huge post-surge increase in aquaculture in AP (vs TN): about +4 percentage points.
# Odisha sees a 0.91 pp smaller increase in aquaculture post-surge than TN.
# higher elevation slightly reduces aquaculture (likely due to lower salinity or reduced inundation risk).

# Before storms, aquaculture is concentrated in high-saline villages. 
# After a storm surge, aquaculture rises more in these already saline areas, suggesting that saline villages may double down on aquaculture rather than revert to traditional agriculture. 
# Low-saline villages show no significant response to storms. This reinforces the idea that salinity creates a path dependency for aquaculture expansion.
# Andhra Pradesh continues to show the largest post-surge increases, hinting at policy support or infrastructure enabling rapid aquaculture expansion.
# Odisha, by contrast, may face constraints or damage that depress post-surge growth. Could be owing to structural poverty. 

# Although, the r2 is less than it is with persistent salinity


# Print results of different saline conditcions 
# All with - FE (Distance) + Clustered SEs + Interaction terms 
models <- list()
models[['Salinity Area']] <- model_3
models[['Persistent Salinity']] <- model_6
models[['Above median Saline']] <- model_9

msummary(models, stars = c('*' = .1, '**' = .05, '***' = .01),gof_omit=c("BIC|AIC|RMSE|R2 Within Adj."),coef_omit=c("(Intercept)"), filename = 'table.rtf')



# Two-way fixed effects 
# Note: Two-way fixed effect controls for village-level fixed effects (absorbing time-invariant differences between villages) - which is why NEAR_DIST and DEM_avg get dropped.
# It controls for year fixed effects (accounting for changes common to all villages in a given year).
# The model focuses on within-village variation over time rather than cross-sectional differences.

data_agg <- data %>%
  group_by(UniqueID, Year) %>%
  summarize(
    Aqua_ha = mean(Aqua_ha, na.rm = TRUE),
    Aqua_perc = mean(Aqua_perc, na.rm = TRUE),
    Saline_perc_norm = mean(Saline_perc_norm, na.rm = TRUE),
    postSurge = max(postSurge, na.rm = TRUE), # A village is affected if any record says so
    salinity_storm_interaction <- Saline_perc_norm * postSurge,
    DEM_avg = mean(DEM_avg, na.rm = TRUE),
    NEAR_DIST = mean(NEAR_DIST, na.rm = TRUE),

    .groups = "drop"
  )

summary(data_agg)

# Set up panel data structure
pdata <- pdata.frame(data_agg, index = c("UniqueID", "Year"))

# Model 10: Fixed effects model (two-way fixed effect - village and year fixed effects)
model_10 <- plm(Aqua_perc ~ postSurge + Saline_perc_norm + 
                 postSurge * Saline_perc_norm, 
               data = pdata, 
               model = "within", 
               effect = "twoways")
summary(model_10)

# After a storm surge, aquaculture increases by ~0.27 pp, holding salinity constant.
# In pre-surge years, a one-unit increase in normalized salinity reduces aquaculture by 0.27 pp.
# After a storm, higher salinity is associated with more aquaculture: the negative effect of salinity flips post-surge.
# => So, after a storm surge, salinity no longer deters aquaculture; in fact, the relationship flips to slightly positive, though small.
# After storm-induced salinization, farmers may switch to aquaculture in saline lands.

# The R2 may seem low, but it is expected because two-way FE models absorb a lot of variation through FE. 
# Besides, this model explain within-village, within-year variation, and not between village or across year variations. 
# Patterns in coefficients being significant is key here. 

# After a storm surge, salinity becomes a trigger for aquaculture — possibly because traditional crops can no longer grow. 
# This supports the idea that post-surge salinization causes land use change toward aquaculture, particularly in areas that become newly saline.


# Testing for serial correlation
pbgtest(model_10)

# This tests whether the error terms within a unit over time (e.g., for a village) are correlated.
# Serial correlation violates a key assumption of standard regression models

# Since serial correlation is present, we could: 
# 1. Use robust standard errors with vcovHC()
# 2. Consider including lagged dependent variables in the model 
# 3. Consider using first differences or other dynamic panel data methods
# 4. Try Newey-West standard errors or other corrections

# Using robust standard errors
robust_model_10 <- coeftest(model_10, vcov = vcovHC(model_10, type = "HC1"))
print(robust_model_10)
# The effect estimates remain mostly the same and highly statistically significant, despite accounting for serial correlation and heteroskedasticity. 


models <- list()
models[['Robust Two-way FE']] <- robust_model_10

msummary(models, stars = c('*' = .1, '**' = .05, '***' = .01),gof_omit=c("BIC|AIC|RMSE|R2 Within Adj."),coef_omit=c("(Intercept)"), filename = 'table.rtf')
# All three effects remain statistically significant at the 0.001 level, even with robust standard errors.
# Before surge: Salinity reduces aquaculture — each unit of Saline_perc_norm is associated with ~0.27 pp lower aquaculture.
# After surge: The interaction offsets this effect — total effect of salinity becomes:
# − 0.266 +  0.275=  +  0.009
# No negative effect anymore — salinity no longer deters aquaculture.

# => Villages that become saline after a storm are more likely to switch to aquaculture.
# The interaction effect is both statistically and substantively important, even after accounting for village and year fixed effects and correcting for serial correlation.






# ----------------------
# VISUALIZATIONS
# ----------------------

# Model 3 visualisation - moving FE to a predictor to enable prediction 
data$State <- relevel(factor(data$State), ref = "TN")

model_3 <- feols(Aqua_perc ~ postSurge + Saline_perc_norm + 
                   salinity_storm_interaction + (State * postSurge) + DEM_avg + Sea_Dist, 
                 data = data, 
                 vcov = ~UniqueID ) 
summary(model_3)

# Visualizing Predicted Values by State and postSurge
# Step 1: Create a new data frame for prediction
# Get means of continuous covariates
mean_salinity <- mean(data$Saline_perc_norm, na.rm = TRUE)
mean_dem <- mean(data$DEM_avg, na.rm = TRUE)

# Create prediction grid
pred_data <- expand.grid(
  postSurge = c(0, 1),
  State = c("AP", "OD", "TN"),
  Saline_perc_norm = mean_salinity,
  DEM_avg = mean_dem,
  Sea_Dist = unique(data$Sea_Dist)[1],  # pick one for fixed effects
  salinity_storm_interaction = pred_data$postSurge * pred_data$Saline_perc_norm
)

# Step 2: Generate predicted values using predict()
# Generate predicted values (with confidence intervals)
pred_data$pred <- predict(model_3, newdata = pred_data, se = TRUE)$fit
pred_data$se <- predict(model_3, newdata = pred_data, se = TRUE)$se.fit

# Confidence intervals
pred_data$lower <- pred_data$pred - 1.96 * pred_data$se
pred_data$upper <- pred_data$pred + 1.96 * pred_data$se


# Step 3: Plot the predicted values
ggplot(pred_data, aes(x = factor(postSurge), y = pred, fill = State)) +
  geom_col(position = position_dodge(0.6), width = 0.5) +
  geom_errorbar(aes(ymin = lower, ymax = upper), position = position_dodge(0.6), width = 0.2) +
  labs(
    x = "Post Surge (0 = Before, 1 = After)",
    y = "Predicted Aquaculture %",
    title = "Predicted Aquaculture Share by State and Storm Surge",
    fill = "State"
  ) +
  theme_minimal()



# Robust Model 10 Visualization 
# Plot predicted Aqua_perc across a range of Saline_perc_norm values using two-way FE (village and year) model. 
# In regression analysis, especially when using models with interactions or nonlinearities, the marginal effect tells us: 
# "How much does the dependent variable change when one independent variable increases by a small amount, holding everything else constant?"
# In this case, since we are using an interaction in the equation: 
# Aqua_perc = b0 + b1 Saline_perc_norm + b2 postSurge + b3 (postSurge * Saline_perc_norm) 
# dAqua_perc/dsaline = b1 + b3 x postSurge (differential eqn)
# then b1 is the baseline sensitivity of aquaculture to salinity and b3 is additional sensitivity due to storm 
# Visualizing this : 

# Step 1: Generate a smooth sequence of salinity values and create prediction data by State and Surge status
sal_seq <- seq(min(data$Saline_perc_norm, na.rm = TRUE),
               max(data$Saline_perc_norm, na.rm = TRUE), 
               length.out = 100)

# Create all combinations of salinity x postSurge x State
states <- unique(data$State)
marg_data <- expand.grid(
  Saline_perc_norm = sal_seq,
  postSurge = c(0, 1),
  State = states
)

# Add interaction term
marg_data$postSurge_Saline <- marg_data$postSurge * marg_data$Saline_perc_norm

# Step 3: Get robust coefficients and variance-covariance matrix
coefs <- coef(robust_model_10) |> as.numeric()
vcov_mat <- vcovHC(model_10, type = "HC1")


# Design matrix
X <- model.matrix(~ postSurge + Saline_perc_norm + postSurge_Saline - 1, data = marg_data)
# The -1 removes the intercept

dim(X)        # (n, k)
length(coefs) # should be k


# Predictions and standard errors
marg_data$predicted <- as.numeric(X %*% coefs)
marg_data$se <- sqrt(rowSums((X %*% vcov_mat) * X))


# Confidence intervals (95%)
marg_data <- marg_data %>%
  mutate(
    lower = predicted - 1.96 * se,
    upper = predicted + 1.96 * se,
    SurgePeriod = ifelse(postSurge == 1, "Post-surge", "Pre-surge")
  )

# Step 4: Plot
ggplot(marg_data, aes(x = Saline_perc_norm, y = predicted, color = SurgePeriod, fill = SurgePeriod)) +
  geom_line(size = 1.2) +
  geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.2, linetype = 0) +
  facet_wrap(~ State) +
  labs(
    title = "Marginal Effect of Salinity on Aquaculture by State",
    x = "Normalized Salinity (%)",
    y = "Predicted Aquaculture (%)",
    color = "Storm Surge Period",
    fill = "Storm Surge Period"
  ) +
  theme_minimal(base_size = 14)

# The downward pre-surge suggests negative marginal effect of salinity 
# the upward post-surge suggests that surges reverse the penalty of salinity on aquaculture 
# The gap between lines increases with salinity → storm surges shift the relationship especially for more saline villages.


# Option 2: plot group means instead of fitted lines 
# Bin salinity for grouping (e.g., 50 bins across full range)
data_binned <- data %>%
  mutate(
    Salinity_bin = cut(Saline_perc_norm, breaks = 10, include.lowest = TRUE),
    SurgePeriod = ifelse(postSurge == 1, "Post-surge", "Pre-surge")
  ) %>%
  group_by(Salinity_bin, postSurge, SurgePeriod, State) %>%
  summarise(
    Saline_perc_norm = mean(Saline_perc_norm, na.rm = TRUE),  # bin center
    mean_aqua = mean(Aqua_perc, na.rm = TRUE),
    se = sd(Aqua_perc, na.rm = TRUE) / sqrt(n()),
    n = n(),
    .groups = "drop"
  ) %>%
  mutate(
    lower = mean_aqua - 1.96 * se,
    upper = mean_aqua + 1.96 * se
  )

# Calculate mean salinity (for vertical line)
sal_mean <- mean(data$Saline_perc_norm, na.rm = TRUE)

# Plot group means with CI ribbons
ggplot(data_binned, aes(x = Saline_perc_norm, y = mean_aqua,
                        color = SurgePeriod, fill = SurgePeriod)) +
  geom_line(size = 1.2) +
  geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.2, linetype = 0) +
  geom_vline(xintercept = sal_mean, linetype = "dashed", color = "black") +
  facet_wrap(~ State) +
  labs(
    title = "Predicted Aquaculture by Salinity, State and Storm Surge",
    x = "Normalized Salinity (%)",
    y = "Predicted Aquaculture (%)",
    color = "Storm Surge Period",
    fill = "Storm Surge Period"
  ) +
  theme_minimal(base_size = 14)






# Compare the first set of regression coefficients 
a<-read.csv("outputs/RQ1_2013-25_AquaHa_Storm_Salinity_Models.csv")

# Use standard errors to calculate CIs
a$upr <- a$est+(1.96*a$se)
a$lwr <- a$est-(1.96*a$se)

# Relevel the first grouping variable so that the plot is ordered the way I want
head(a)
a$Variable <- as.factor(a$Variable)
levels(a$Variable)
a$Variable <- factor(a$Variable, levels=c("Salinity (5 yrs)","Saline Area", "Surge"))
a$Controls<-as.factor(a$Controls)


# Restrict the data only to estimates from select models
a1 <- subset(a, Controls == "OLS + Interaction" | Controls == "FE (State + Sea Distance)" | Controls == "Persistent Salinity (OLS + Interaction)" | Controls == "Persistent Salinity (FE)")

# Make sure that the Controls column is a factor and reorder as needed
levels(a1$Controls)
a1$Controls <- factor(a1$Controls, levels=c("OLS + Interaction","FE (State + Sea Distance)","Persistent Salinity (OLS + Interaction)", "Persistent Salinity (FE)"))
a1$Controls<-as.factor(a1$Controls)

# Run the plot
# This is the plot with the horizontal lines
plot<-ggplot(a1,aes(est,Variable,group=Controls))+
  geom_errorbarh(aes(xmax = upr, xmin = lwr),height=0.5,position=position_dodge(width=0.5), colour="black")+
  geom_point(size=5, aes(shape=Controls,color=Controls),position=position_dodge(width=0.5))+
  theme_bw()+
  ggtitle("Relationship between Aquaculture, Surge & Salinity across model specs")+
  ylab("")+   theme(legend.position="right")+
  theme(legend.text = element_text(size = 12),legend.title = element_text(size = 12))+
  theme(plot.title = element_text(hjust = 0.5),panel.border = element_rect(colour = "black",linewidth=1))+
  theme(axis.title = element_text(size = 16,face="bold"),axis.text=element_text(size=16,color="black"))+
  theme(plot.title = element_text(size = 16,face = "bold"))+xlab("") +
  theme(panel.grid.major = element_line(colour="grey", linewidth=0.5)) +
  theme(panel.grid.minor = element_line(colour="grey", linewidth=0.5))+
  theme(panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        panel.border = element_blank(),
        axis.line = element_line(color = 'black'))+theme(legend.title = element_blank(), 
                                                         legend.key.size = unit(3, 'cm'))+
  geom_hline(yintercept = 4.5, linetype = "dotted",color ="black",linewidth = 1)+
  geom_hline(yintercept = 6.5, linetype = "dotted",color ="black",linewidth = 1)+
  geom_vline(xintercept = 0, linetype = "longdash",size = 0.5)+
  
  scale_color_manual(values = c( "#ffab33","#990000", "#666633", "#003300"))+
  scale_shape_manual(values = c(16,15,18,17))

plot(plot)

ggsave("outputs/RQ1_2013-25_AquaHa_Storm_Salinity_Models.png",plot = plot, width = 12, height = 6, units = "in",dpi=300)





# Visualize the interaction effect
# Create predicted values based on Model 2
grid_data <- expand.grid(
  Saline_perc_norm = seq(min(data$Saline_perc_norm, na.rm = TRUE),
                            max(data$Saline_perc_norm, na.rm = TRUE),
                            length.out = 100),
  postSurge = c(0, 1),
  NEAR_DIST = mean(data$NEAR_DIST, na.rm = TRUE),
  DEM_avg= mean(data$DEM_avg, na.rm = TRUE)
)

grid_data$salinity_storm_interaction <- grid_data$Saline_perc_norm * grid_data$postSurge
grid_data$predicted <- predict(model_2, newdata = grid_data)

# Plot the interaction effect
ggplot(grid_data, aes(x = Saline_perc_norm, y = predicted, color = factor(postSurge))) +
  geom_line() +
  labs(title = "Interaction Effect of Salinity and Storm Surge on Aquaculture",
       x = "Saline Area Percentage",
       y = "Predicted Aquaculture",
       color = "Storm Surge") +
  theme_minimal()



# Map Flood - Saline areas 
# map1 <- merge(map_dist_shp, dist_state_data, by="District_name", all.x=T)
# 
# # Poverty 
# summary(map1$poverty_share) 
# # Some districts are outliers i.e. have too high poverty (>0.8) may want to remove those for mapping? 
# map_poverty <- ggplot() +
#   geom_sf(data=map_dist_shp, colour = "lightgrey", fill="lightgrey")+ # To map all districts even if data not available
#   geom_sf(data=map1, aes(fill = poverty_share), colour=alpha("dark grey", 1/2), size=0.2) +
#   geom_sf(data=map_state_shp, colour="white", fill=NA, size = 4) + # draw polygon outline of states and fill nothing
#   coord_sf() +
#   theme_minimal() +
#   ggtitle("Poverty share across districts: % of people in 0-20% MPCE Quantile") +
#   theme(axis.line = element_blank(), axis.text = element_blank(),
#         axis.ticks = element_blank(), axis.title = element_blank(), 
#         plot.title = element_text(size = 12),
#         legend.key.size = unit(2, "lines"),
#         panel.grid.major = element_blank(),
#         panel.grid.minor = element_blank()) +
#   scale_fill_viridis(option="cividis", breaks = seq(0, 1, by = 0.2))
# # Adjust legend position and size
# map_poverty <- map_poverty +
#   theme(legend.position = "bottom",  # Adjust legend position
#         legend.text = element_text(size = 10),  # Adjust legend text size
#         legend.title = element_text(size = 10))  # Adjust legend title size
# #Export map
# ggsave("Map_poverty.png",plot = map_poverty, width = 8, height = 8, units = "in",dpi=600)
# 
