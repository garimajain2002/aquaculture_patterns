
# RQ: How does aquaculture expansion trend vary after a place is affected by a surge? 
# - As compared to its previous growth (within-place event study)
# - As compared to the counter factual (matched dif-in-dif) 


# ----------------------
# Data and prep 
# ----------------------
library(fixest)
library(tidyr)    
library(tidyverse)
library(zoo) #for rollmean() for calculating smoother percentage changes over 3-5 years
library(slider)
library(modelsummary)
library(ggplot2)
library(dplyr)
library(did)


getwd()

df <- read.csv(unz("data/aqua_salinity_surge_1990-2025.zip", 
                     "aqua_salinity_surge_1990-2025.csv"))

# Make TN reference state
df$State <- relevel(factor(df$State), ref = "TN")

# Create dummy for treated villages (those ever affected)
df <- df %>%
  group_by(UniqueID) %>%
  mutate(treated = any(postSurge == 1)) %>%
  ungroup()

# Add a flag variable for any key missing variables 
# If any of these is NA, make flag = NA : Aqua_perc , postSurge , avg_salinity_5yr , Lag_Aqua
df$flag <- ifelse(rowSums(is.na(df[c("Aqua_perc", "postSurge", "avg_salinity_5yr", "Lag_Aqua")])) > 0, NA, 0)

table(df$flag, useNA = "always")  # Shows count of flagged vs non-flagged rows
sum(is.na(df$flag))   # Total number of rows with missing variables


# --------------------------------
# Dif-in-Dif model (All+By States)
# --------------------------------

did_model <- feols(
  Aqua_perc ~ avg_salinity_5yr + i(postSurge, treated, ref = 0) + flag | UniqueID + Year,
  data = df
)
summary(did_model)

models <- list()
models[['Diff-in-Diff']] <- did_model
msummary(models, stars = c('*' = .1, '**' = .05, '***' = .01),gof_omit=c("BIC|AIC|RMSE|R2 Within Adj."),coef_omit=c("(Intercept)"), filename = 'table.rtf')

# Interpretation: 
# Each 1% increase in avg. salinity over past 5 years is associated with a 0.201 pp increase in aquaculture share. 
# Villages affected by surge show a +0.244 increase in aquaculture after the surge compared to unaffected villages. This is statistically significant at 99% CI (p < 0.01)

# This supports a causal interpretation: storm surges as well as salinity push aquaculture expansion.
# The effect (+0.148 pp) may seem modest, but important if aggregated across space and time. 


# DiD with lag_aqua also (reflecting the same set as the main panel regression models)
did_model_2 <- feols(
  Aqua_perc ~ i(postSurge, treated, ref = 0) + avg_salinity_5yr + Lag_Aqua + flag | UniqueID + Year,
  data = df
)
summary(did_model_2)

models <- list()
models[['Diff-in-Diff']] <- did_model_2
msummary(models, stars = c('*' = .1, '**' = .05, '***' = .01),gof_omit=c("BIC|AIC|RMSE|R2 Within Adj."),coef_omit=c("(Intercept)"), filename = 'table.rtf')

# Interpretation: 
# Villages affected by surge show a +0.121 increase in aquaculture after the surge compared to unaffected villages.
# Each 1% increase in avg. salinity over past 5 years is associated with a 0.068 pp increase in aquaculture share. 
# These results are about the same as the regular regression model from earlier.  



# DiD model for States 
# Pooled model with state interactions
pooled_model <- feols(
  Aqua_perc ~ 
    avg_salinity_5yr * State + 
    i(postSurge, treated, ref = 0) * State + flag | 
    UniqueID + Year,
  cluster = ~UniqueID,
  data = df
)

summary(pooled_model)

models <- list()
models[['Pooled Diff-in-Diff']] <- pooled_model
msummary(models, stars = c('*' = .1, '**' = .05, '***' = .01),gof_omit=c("BIC|AIC|RMSE|R2 Within Adj."),coef_omit=c("(Intercept)"), filename = 'table.rtf')

# Effect of salinity : TN (baseline): -0.13 (salinity reduces aquaculture); AP: -0.13 + 0.62 = +0.49; OD: -0.13 + 0.29 = +0.16
# Effect of surges: TN: Full effect = +3.35 pp (strong, positive, significant); AP: +1.02 pp; OD: Effect dropped due to collinearity — this likely means there's no variation left (e.g., no untreated control units or no surge timing variation within OD) once FE and other interactions are absorbed.
# + State indicators are collinear with fixed effects 

# Check variation in surge timing within OD 
od_df <- df %>% filter(State == "OD")
table(od_df$treated, useNA = "ifany")

od_df %>%
  filter(treated == 1) %>%
  distinct(UniqueID, first_surge_year) %>%
  count(first_surge_year)

od_df %>%
  group_by(Year, treated) %>%
  summarize(
    n = n(),
    avg_post = mean(postSurge, na.rm = TRUE)
  ) %>%
  filter(treated == 1)

# Majority of villages in OD were treated in two surge years - 1971 and 1999.
# When 1971, most treated villages are already in post-surge period 
# The remaining (1999 surge) account for the rest, flipping postSurge to 1 only in 1999 onward.
# No variation over time within a village, and little variation between villages, especially within year.



# -------------------------------------------------------------------------------------------------------
# Dynamic Dif-in-Dif model with staggard treatment timing (Callaway & Sant'Anna / Sun and Abraham method)
# -------------------------------------------------------------------------------------------------------

# 1. Drop villages treated before 1990 from the analysis (satisfying the assumption that no one is treated in t=1)
# 2. All untreated villages, G = inf.
# 3. Since states differ substantially in trends, treatment timing and treatment effect heterogeneity, running separate DiD models for each state can improve validity and interpretability. or include State FE with State x Time interaction? Try both. 
# 4. Covariates for matching: 
# Saline (binary for above or below median salinity) OR Saline_perc_norm as dbl variable with actual normalised values, 
# Sea_Dist (categorical variable for distance from the sea) OR NEAR_DIST (cont. dbl variable for actual distance from the sea), 
# DEM_avg (average elevation from the sea), 
# agriculture_percent (% area under agriculture), 
# urban_percent (% area urban), 
# TOT_WORK_P (total working population which can be converted to a % by dividing with TOT_P or total population), and 
# MARGWORK_P (total marginal workers and can be turned into a % by dividing with TOT_WORK_P as a %marginal workers of total working population).



# Define event time for all places (treatment and control)
df$YearSurge = df$first_surge_year

table(df$YearSurge)

# Define G (group) based on treatment time and inf when untreated
df <- df %>%
  mutate(
    G = if_else(is.na(YearSurge),       # those never hit
                Inf,                    # → G = ∞
                YearSurge)              # otherwise G = first surge year
  )


table(df$G)

# Keep only first treated in >=1990. Drop the rest. 
df_1990 <- df %>%
  filter(G > 1990 | is.infinite(G))

# Make UniqueID numeric
df_1990 <- df_1990 %>%
  mutate(numeric_id = as.numeric(factor(UniqueID)))

# Check for NAs in covariates before using 
df_1990 %>%
  select(Aqua_perc, Year, UniqueID, G, flag,
         NEAR_DIST, Saline_perc_norm, DEM_avg,
         agriculture_percent, urban_percent,
         TOT_WORK_P, TOT_P, MARGWORK_P) %>%
  summarise_all(~ sum(is.na(.)))

# Compute some covariates
df_1990 <- df_1990 %>%
  mutate(
    WORK_PERC = TOT_WORK_P / TOT_P,
    MARG_WORK_PERC = MARGWORK_P / TOT_WORK_P
  )



# One state at a time
# OD 
df_OD <- df_1990 %>% filter(State == "OD")

table(df_OD$G)


# Combined DiD using Callaway and Sant'Anna (2020) did package 
out_OD <- att_gt(
  yname    = "Aqua_perc",
  gname    = "G",
  idname   = "numeric_id",
  tname    = "Year",
  xformla  = ~ NEAR_DIST + Saline_perc_norm + DEM_avg +
    agriculture_percent + urban_percent +
    WORK_PERC + MARG_WORK_PERC,
  data     = df_OD,
  control_group = "nevertreated",     # units with G = Inf
  est_method    = "dr",                  # doubly‐robust
  panel = FALSE                    # FALSE to allow unbalanced panel - required due to dropped early-treated units
)

# Event‐study plot
ggdid(out_OD)

p <- ggdid(out_OD)

ggsave("visuals/DiD/OD_event_study.jpg", plot = p, width = 14, height = 6, dpi = 300)


# ATT estimates 
summary(out_OD)

# Average treatment effect on the treated
agg_OD <- aggte(out_OD, type = "group", na.rm = TRUE) #Aggregate Group-Time Average Treatment Effects
summary(agg_OD)



# Alternate plot for OD

# Step 1: Tidy the object
att_df_clean <- tidy(out_OD)

# Step 2: Create event time and classify pre/post
att_df_clean$event_time <- att_df_clean$time - att_df_clean$group
att_df_clean$period <- ifelse(att_df_clean$event_time < 0, "Pre", "Post")

# Step 3: Plot
p<- ggplot(att_df_clean, aes(x = event_time, y = estimate, color = period)) +
  geom_point() +
  geom_line(aes(group = 1)) +
  geom_errorbar(aes(ymin = estimate - 1.96 * std.error,
                    ymax = estimate + 1.96 * std.error),
                width = 0.2) +
  facet_wrap(~ group, scales = "free_x") +
  scale_color_manual(values = c("Pre" = "orange", "Post" = "steelblue")) +
  labs(title = "Event Study: Pre vs Post Treatment ATT",
       x = "Event Time (Year - Group Year)", y = "ATT",
       color = "Period") +
  theme_minimal()

plot(p) 

ggsave("visuals/DiD/OD_event_study_alt.jpg", plot = p, width = 14, height = 6, dpi = 300)





# AP
df_AP <- df_1990 %>% filter(State == "AP")

table(df_AP$G)

# COMBINED TAKES A LONG TIME
# Combined DiD using Callaway and Sant'Anna (2020) did package 
out_AP <- att_gt(
  yname    = "Aqua_perc",
  gname    = "G",
  idname   = "numeric_id",
  tname    = "Year",
  xformla  = ~ NEAR_DIST + Saline_perc_norm + DEM_avg +
    agriculture_percent + urban_percent +
    WORK_PERC + MARG_WORK_PERC,
  data     = df_AP,
  control_group = "nevertreated",     # units with G = Inf
  est_method    = "dr",                  # doubly‐robust
  panel = FALSE                    # FALSE to allow unbalanced panel - required due to dropped early-treated units
)

# Event‐study plot
ggdid(out_AP)

p <- ggdid(out_AP)

ggsave("visuals/DiD/AP_event_study.jpg", plot = p, width = 14, height = 6, dpi = 300)



# ATT estimates 
summary(out_AP)

# Average treatment effect on the treated
agg_AP <- aggte(out_AP, type = "group", na.rm = TRUE)
summary(agg_AP)


# Alternate plot for AP

# Step 1: Tidy the object
att_df_clean <- tidy(out_AP)

# Step 2: Create event time and classify pre/post
att_df_clean$event_time <- att_df_clean$time - att_df_clean$group
att_df_clean$period <- ifelse(att_df_clean$event_time < 0, "Pre", "Post")

# Step 3: Plot
p<- ggplot(att_df_clean, aes(x = event_time, y = estimate, color = period)) +
  geom_point() +
  geom_line(aes(group = 1)) +
  geom_errorbar(aes(ymin = estimate - 1.96 * std.error,
                    ymax = estimate + 1.96 * std.error),
                width = 0.2) +
  facet_wrap(~ group, scales = "free_x") +
  scale_color_manual(values = c("Pre" = "orange", "Post" = "steelblue")) +
  labs(title = "Event Study: Pre vs Post Treatment ATT",
       x = "Event Time (Year - Group Year)", y = "ATT",
       color = "Period") +
  theme_minimal()

plot(p) 

ggsave("visuals/DiD/AP_event_study_alt.jpg", plot = p, width = 14, height = 6, dpi = 300)






# TN
df_TN <- df_1990 %>% filter(State == "TN")

table(df_TN$G)
table(df_TN$Year, df_TN$G)  # To see coverage per year per group

# Combined DiD using Callaway and Sant'Anna (2020) did package 
out_TN <- att_gt(
  yname    = "Aqua_perc",
  gname    = "G",
  idname   = "numeric_id",
  tname    = "Year",
  xformla  = ~ NEAR_DIST + Saline_perc_norm + DEM_avg +
    agriculture_percent + urban_percent +
    WORK_PERC + MARG_WORK_PERC,
  data     = df_TN,
  control_group = "nevertreated",     # units with G = Inf
  est_method    = "dr",                  # doubly‐robust
  panel = FALSE                    # FALSE to allow unbalanced panel - required due to dropped early-treated units
)


# Event‐study plot
ggdid(out_TN)

p <- ggdid(out_TN)

ggsave("visuals/DiD/TN_event_study.jpg", plot = p, width = 14, height = 18, dpi = 300)

# ATT estimates 
summary(out_TN)

# Average treatment effect on the treated
agg_TN <- aggte(out_TN, type = "group", na.rm = TRUE)
summary(agg_TN)




# Alternate plot for TN without 1992 

# Remove the group == 1992 since all NA
# Step 1: Tidy the object and remove 1992
att_df <- tidy(out_TN)
att_df_clean <- att_df[att_df$group != 1992, ]

# Step 2: Create event time and classify pre/post
att_df_clean$event_time <- att_df_clean$time - att_df_clean$group
att_df_clean$period <- ifelse(att_df_clean$event_time < 0, "Pre", "Post")

# Step 3: Plot
p<- ggplot(att_df_clean, aes(x = event_time, y = estimate, color = period)) +
  geom_point() +
  geom_line(aes(group = 1)) +
  geom_errorbar(aes(ymin = estimate - 1.96 * std.error,
                    ymax = estimate + 1.96 * std.error),
                width = 0.2) +
  facet_wrap(~ group, scales = "free_x") +
  scale_color_manual(values = c("Pre" = "orange", "Post" = "steelblue")) +
  labs(title = "Event Study: Pre vs Post Treatment ATT",
       x = "Event Time (Year - Group Year)", y = "ATT",
       color = "Period") +
  theme_minimal()

plot(p) 

ggsave("visuals/DiD/TN_event_study_alt.jpg", plot = p, width = 14, height = 6, dpi = 300)



# Combine all event studies into one graphic (alternate visualization) 
# Tidy and label each output
tidy_TN  <- tidy(out_TN)  %>% filter(group != 1992) %>% mutate(state = "Tamil Nadu")
tidy_AP  <- tidy(out_AP)  %>% mutate(state = "Andhra Pradesh")
tidy_OD  <- tidy(out_OD)  %>% mutate(state = "Odisha")

# Combine all
att_combined <- bind_rows(tidy_TN, tidy_AP, tidy_OD)

# Create event time and pre/post indicator
att_combined <- att_combined %>%
  mutate(
    event_time = time - group,
    period = ifelse(event_time < 0, "Pre", "Post"),
    group_label = paste(state, group)
  )

# Plot with facets by state and cohort (group year)
p_combined <- ggplot(att_combined, aes(x = event_time, y = estimate, color = period)) +
  geom_point(size = 2) +
  geom_line(aes(group = 1)) +
  geom_errorbar(aes(ymin = estimate - 1.96 * std.error,
                    ymax = estimate + 1.96 * std.error),
                width = 0.3) +
  facet_wrap(~ group_label, ncol = 2, scales = "free_x") +
  scale_color_manual(values = c("Pre" = "orange", "Post" = "steelblue")) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  labs(
    title = "Dynamic Effects of Storm Surge Exposure on Aquaculture by Event Cohort",
    x = "Event Time (Years Since Surge)", y = "ATT (pp effect)",
    color = "Period"
  ) +
  theme_minimal(base_size = 14)

# Display
plot(p_combined)

# Save
ggsave("visuals/DiD/Combined_EventStudy_ATT.jpg", plot = p_combined, width = 12, height = 8, dpi = 300)







# Compute overall ATT using OD, AP, TN ATTs
# Extract group-level ATT and SE for each state
att_od  <- agg_OD$overall.att
se_od   <- agg_OD$overall.se
n_od <- df_1990 %>%
  filter(State == "OD") %>%
  filter(is.finite(G)) %>%        # only treated units (G ≠ Inf)
  distinct(UniqueID) %>%
  nrow()

att_ap  <- agg_AP$overall.att
se_ap   <- agg_AP$overall.se
n_ap <- df_1990 %>%      # treated units in AP
  filter(State == "AP") %>%
  filter(is.finite(G)) %>%        # only treated units (G ≠ Inf)
  distinct(UniqueID) %>%
  nrow()  

att_tn  <- agg_TN$overall.att
se_tn   <- agg_TN$overall.se
n_tn <- df_1990 %>%     # treated units in TN
  filter(State == "TN") %>%
  filter(is.finite(G)) %>%
  distinct(UniqueID) %>%
  nrow()  

# Combine into vectors
atts <- c(att_od, att_ap, att_tn)
weights <- c(n_od, n_ap, n_tn)

# Weighted average ATT
overall_att <- weighted.mean(atts, weights)
overall_att


# Combine SEs into variances
vars <- c(se_od^2, se_ap^2, se_tn^2)

# Weighted standard error of the mean
overall_se <- sqrt( sum((weights^2 * vars)) / (sum(weights)^2) )
overall_se

# 95% Confidence Interval
ci_lower <- overall_att - 1.96 * overall_se
ci_upper <- overall_att + 1.96 * overall_se
c(ci_lower, ci_upper)

# Overall plot option 
# Get tidy group-time estimates from each state
tidy_OD <- tidy(out_OD) %>% mutate(State = "OD")
tidy_AP <- tidy(out_AP) %>% mutate(State = "AP")
tidy_TN <- tidy(out_TN) %>% mutate(State = "TN")

# Combine into a single dataframe
att_all <- bind_rows(tidy_OD, tidy_AP, tidy_TN)

# Remove 1992 since NA
att_all_clean <- att_all %>%
  filter(!is.na(estimate)) %>%
  filter(group != 1992)  # optional, based on your note earlier

# Plot
ggplot(att_all_clean, aes(x = time, y = estimate, color = State)) +
  geom_point() +
  geom_line() +
  geom_errorbar(aes(ymin = estimate - 1.96 * std.error,
                    ymax = estimate + 1.96 * std.error),
                width = 0.2) +
  facet_wrap(~ group, scales = "free_x") +  # separate panel per event group
  geom_hline(yintercept = 0, linetype = "dashed") +
  labs(title = "Group-Time ATT Estimates by State",
       x = "Year", y = "ATT", color = "State") +
  theme_minimal(base_size = 14)

ggsave("visuals/DiD/All_event_study.jpg", width = 14, height = 6, dpi = 300)









# Optional : Sun & Abraham (2021) single dynamic event-study plot 
# -----

# Group Size 
n_OD <- sum(!is.infinite(out_OD$G))
n_AP <- sum(!is.infinite(out_AP$G))
n_TN <- sum(!is.infinite(out_TN$G))

# Extract dynamic event-study estimates from each
# OD
agg_OD <- aggte(out_OD, type = "dynamic", na.rm = TRUE)
dyn_OD <- data.frame(
  event_time = agg_OD$egt,
  estimate = agg_OD$att.egt,
  std.error = agg_OD$se.egt,
  State = "OD",
  N = n_OD
)


# AP
agg_AP <- aggte(out_AP, type = "dynamic", na.rm = TRUE)
dyn_AP <- data.frame(
  event_time = agg_AP$egt,
  estimate = agg_AP$att.egt,
  std.error = agg_AP$se.egt,
  State = "AP",
  N = n_AP
)


# TN
agg_TN <- aggte(out_TN, type = "dynamic", na.rm = TRUE)
dyn_TN <- data.frame(
  event_time = agg_TN$egt,
  estimate = agg_TN$att.egt,
  std.error = agg_TN$se.egt,
  State = "TN",
  N = n_TN
)

# Combine dataframes 
dyn_all <- bind_rows(dyn_OD, dyn_AP, dyn_TN)


# Average across states for each event time
dyn_avg <- dyn_all %>%
  group_by(event_time) %>%
  summarize(
    estimate = mean(estimate, na.rm = TRUE),
    std.error = sqrt(mean(std.error^2, na.rm = TRUE) / 3),  # average of squared SEs
    .groups = "drop"
  )


ggplot(dyn_avg, aes(x = event_time, y = estimate)) +
  geom_line(color = "#0072B2", size = 1.2) +
  geom_point(color = "#0072B2", size = 2.5) +
  geom_errorbar(aes(ymin = estimate - 1.96 * std.error,
                    ymax = estimate + 1.96 * std.error),
                width = 0.3, color = "#0072B2") +
  geom_vline(xintercept = 0, linetype = "dashed", color = "black") +
  geom_hline(yintercept = 0, linetype = "dotted", color = "grey40") +
  labs(
    title = "Average Dynamic ATT across All Storm Events and States",
    x = "Years Relative to First Storm Surge Event",
    y = "Average ATT on Aquaculture (%)"
  ) +
  theme_minimal(base_size = 14)

ggsave("outputs/Overall_ATT_SunAbraham.jpg", width = 14, height = 6, dpi = 300)






# TAKES A LONG TIME TO RUN 
# ----
# Compute all together to get the overall ATT 
out_all <- att_gt(
  yname    = "Aqua_perc",
  gname    = "G",
  idname   = "numeric_id",
  tname    = "Year",
  xformla  = ~ NEAR_DIST + Saline_perc_norm + DEM_avg +
    agriculture_percent + urban_percent +
    WORK_PERC + MARG_WORK_PERC,
  data     = df_1990,
  control_group = "nevertreated",
  est_method    = "dr",
  panel = FALSE
)

# Aggregate to get overall ATT
agg_overall <- aggte(out_all, type = "group")
summary(agg_overall)


# Event‐study plot
ggdid(out_all)

p <- ggdid(out_all)

ggsave("outputs/All_event_study.jpg", plot = p, width = 14, height = 18, dpi = 300)

# ATT estimates 
summary(out_all)

# Average treatment effect on the treated
agg_all <- aggte(out_all, type = "group", na.rm = TRUE)
summary(agg_all)







# ------ OLD Script ------

#----------------
# Visualise DiD
#----------------
# Define event time for all places (treatment and control)
df$YearSurge = df$first_surge_year

# Set 1990 for AP, 1999 for OD, and 2004 for TN
df$YearSurge[df$treated == FALSE & df$State == "AP"] <- 1990
df$YearSurge[df$treated == FALSE & df$State == "OD"] <- 1999
df$YearSurge[df$treated == FALSE & df$State == "TN"] <- 1992 # Use 2004 only once tsunami is added in the df

table(df$first_surge_year)
table(df$YearSurge)

df$event_time = df$Year - df$YearSurge

table(df$event_time)

# ! Check whether to work with Aqua_perc vs Aqua_ha


# # Assign some psedo year to control villages, else will get filtered out. 
# # Get treated villages' actual surge years
# treated_years <- df %>%
#   filter(treated == TRUE) %>%
#   distinct(UniqueID, first_surge_year)
# 
# # Sample placebo years for controls
# set.seed(123)
# control_placebo_years <- df %>%
#   filter(treated == FALSE) %>%
#   distinct(UniqueID) %>%
#   mutate(
#     placebo_surge = sample(treated_years$first_surge_year, size = n(), replace = TRUE)
#   )
# 
# # Join back to main df
# df <- df %>%
#   left_join(control_placebo_years, by = "UniqueID") %>%
#   mutate(
#     YearSurge = ifelse(treated == FALSE, placebo_surge, first_surge_year),
#     event_time = Year - YearSurge
#   )

# # Using only interquantile rage
# q5 <- quantile(df$Aqua_perc, 0.05, na.rm = TRUE)
# q95 <- quantile(df$Aqua_perc, 0.95, na.rm = TRUE)

# did_event <- df %>%
#   filter(!is.na(event_time) & event_time >= -12 & event_time <= 12) %>%
#   filter(Aqua_perc >= q5 & Aqua_perc <= q95) %>%
#   group_by(event_time, treated) %>%
#   summarise(mean_aqua = mean(Aqua_perc, na.rm = TRUE), .groups = "drop") %>%
#   mutate(treated = ifelse(treated == 1, "Treated", "Control"))

did_event <- df %>%
  filter(!is.na(event_time) & event_time >= -12 & event_time <= 12) %>%  # e.g., +/- 12 years around treatment
  group_by(event_time, treated) %>%
  summarise(mean_aqua = mean(Aqua_perc, na.rm = TRUE), .groups = "drop") %>%
  mutate(treated = ifelse(treated == 1, "Treated", "Control"))

ggplot(did_event, aes(x = event_time, y = mean_aqua, color = treated)) +
  geom_line(size = 1.2) +
  geom_vline(xintercept = 0, linetype = "dashed") +
  labs(
    title = "Event-Time Plot: Aquaculture Trends Around Storm Surge",
    subtitle = "Year 0 = Year of First Surge per Village",
    x = "Years Since First Storm Surge",
    y = "Mean Aquaculture %",
    color = "Group"
  ) +
  theme_minimal(base_size = 14)



# Test what places have the sudden hike 2-3 years before the event
# Step 1: Identify spike candidate UniqueIDs
spike_ids <- df %>%
  filter(treated == TRUE,
         event_time %in% c(-3, -2, -1)) %>%
  group_by(UniqueID) %>%
  summarise(mean_aqua = mean(Aqua_perc, na.rm = TRUE), .groups = "drop") %>%
  filter(mean_aqua > quantile(df$Aqua_perc, 0.85, na.rm = TRUE)) %>%
  pull(UniqueID)

# Step 2: Return full info for those spike candidates in relevant years
spike_candidates <- df %>%
  filter(UniqueID %in% spike_ids,
         event_time %in% c(-3, -2, -1),
         treated == TRUE) %>%
  select(UniqueID, State, Year, YearSurge, event_time, treated, Aqua_perc)

spike_history <- df %>%
  filter(UniqueID %in% spike_candidates$UniqueID) %>%
  arrange(UniqueID, Year)

ggplot(spike_history, aes(x = event_time, y = Aqua_perc, group = UniqueID)) +
  geom_line(aes(color = State), alpha = 0.3) +
  #geom_smooth(se = FALSE, color = "black") +
  labs(
    title = "Aquaculture Trends: Villages with Pre-Surge Spikes",
    x = "Event Time (Years Relative to Surge)",
    y = "Aquaculture (%)"
  ) +
  theme_minimal() +
  theme(legend.title = element_blank())



# Plot by state 
did_event <- df %>%
  filter(!is.na(event_time) & event_time >= -12 & event_time <= 12) %>%
  group_by(event_time, treated, State) %>%
  summarise(mean_aqua = mean(Aqua_perc, na.rm = TRUE), .groups = "drop") %>%
  mutate(treated = ifelse(treated == 1, "Treated", "Control"))

ggplot(did_event, aes(x = event_time, y = mean_aqua, color = treated)) +
  geom_line(size = 1.2) +
  geom_vline(xintercept = 0, linetype = "dashed") +
  labs(
    title = "Event-Time Plot: Aquaculture Trends Around Storm Surge by State",
    subtitle = "Year 0 = Year of First Surge per Village",
    x = "Years Since First Storm Surge",
    y = "Mean Aquaculture %",
    color = "Group"
  ) +
  facet_wrap(~State) +  # This adds separate panels for each State
  theme_minimal(base_size = 14)




# Check lagged effects (Sun & Abraham (2020))
model_sunab <- feols(
  Aqua_perc ~ avg_salinity_5yr + sunab(YearSurge, Year, ref.p = -1) |
    UniqueID + Year,
  cluster = ~UniqueID,
  data = df
)

iplot(model_sunab, 
      ref.line = -1,
      xlim = c(-10, 10),  # restrict x-axis to ±10
      ci.level = 0.95,
      main = "Event Study: Effect of Surge on Aquaculture",
      xlab = "Years since first surge")


