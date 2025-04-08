# Explore the data and descriptive statistics 
library(dplyr)
library(ggplot2)
library(tidyverse)
library(corrplot)

getwd()

df <- read.csv("outputs/aqua_salinity_surge_2013-2025.csv")
summary(df)

# df of near and very near villages only 
table(df$Sea_Dist)
df_near <- df %>%
  filter(Sea_Dist %in% c("Very Near", "Near"))

summary(df_near)

# df of very near villages only 
table(df$Sea_Dist)
df_vnear <- df %>%
  filter(Sea_Dist %in% c("Very Near"))

summary(df_vnear)



# ################ SUMMMARY STATS ##############

# 1. Mean, Median, Min, Max by year
summary_stats <- df %>%
  group_by(Year) %>%
  summarise(
    mean_aqua = mean(Aqua_perc, na.rm = TRUE),
    median_aqua = median(Aqua_perc, na.rm = TRUE),
    min_aqua = min(Aqua_perc, na.rm = TRUE),
    max_aqua = max(Aqua_perc, na.rm = TRUE),
    mean_saline = mean(Saline_perc, na.rm = TRUE),
    median_saline = median(Saline_perc, na.rm = TRUE),
    min_saline = min(Saline_perc, na.rm = TRUE),
    max_saline = max(Saline_perc, na.rm = TRUE)
  )

print(summary_stats)
# Finding: Overall trend shows increasing aquaculture and salinity over the decade

# Basic Summary statistics by year
yearly_summary <- df %>%
  group_by(Year) %>%
  summarize(
    mean_aquaculture = mean(Aqua_perc, na.rm = TRUE),
    mean_salinity = mean(Saline_perc, na.rm = TRUE),
    storm_affected_villages = sum(postSurge, na.rm = TRUE),
    n_villages = n_distinct(UniqueID)
  )
print(yearly_summary)


# Visualize aquaculture and salinity trends over time
# ALL VILLAGES 
ggplot(yearly_summary, aes(x = Year, y = mean_aquaculture)) +
  geom_line() +
  geom_point() +
  labs(title = "Average Aquaculture Percentage by Year",
       x = "Year", y = "Average Aquaculture (%)") +
  scale_x_continuous(breaks = seq(min(yearly_summary$Year), max(yearly_summary$Year), by = 1)) +  # Ensure integer year labels
  
  theme_minimal()
ggsave("outputs/AquaTrends_2013-25_allVillages.png", width = 8, height = 6, dpi = 300)


ggplot(yearly_summary, aes(x = Year, y = mean_salinity)) +
  geom_line() +
  geom_point() +
  labs(title = "Average Saline Area Percentage by Year",
       x = "Year", y = "Average Saline Area (%)") +
  scale_x_continuous(breaks = seq(min(yearly_summary$Year), max(yearly_summary$Year), by = 1)) +  # Ensure integer year labels
  
  theme_minimal()
ggsave("outputs/SalineTrends_2013-25_allVillages.png", width = 8, height = 6, dpi = 300)



# 2. Mean Aquaculture and Salinity by Storm Surge Status
yearly_summary <- df %>%
  group_by(Year) %>%
  summarize(
    # Mean aquaculture % for affected and not affected villages
    mean_aquaculture_affected = mean(Aqua_perc[postSurge == 1], na.rm = TRUE),
    mean_aquaculture_not_affected = mean(Aqua_perc[postSurge == 0], na.rm = TRUE),
    
    # Standard deviation for aquaculture
    sd_aquaculture_affected = sd(Aqua_perc[postSurge == 1], na.rm = TRUE),
    sd_aquaculture_not_affected = sd(Aqua_perc[postSurge == 0], na.rm = TRUE),
    
    # Mean salinity % for affected and not affected villages
    mean_salinity_affected = mean(Saline_perc[postSurge == 1], na.rm = TRUE),
    mean_salinity_not_affected = mean(Saline_perc[postSurge == 0], na.rm = TRUE),
    
    # Standard deviation for salinity
    sd_salinity_affected = sd(Saline_perc[postSurge == 1], na.rm = TRUE),
    sd_salinity_not_affected = sd(Saline_perc[postSurge == 0], na.rm = TRUE),
    
    # Count of affected and not affected villages
    n_affected = sum(postSurge == 1, na.rm = TRUE),
    n_not_affected = sum(postSurge == 0, na.rm = TRUE),
    
    # Count of storm-affected villages
    storm_affected_villages = sum(postSurge, na.rm = TRUE),
    
    # Total number of unique villages
    n_villages = n_distinct(UniqueID)
  ) %>%
  mutate(
    # Compute standard error (SE) for aquaculture
    se_aquaculture_affected = sd_aquaculture_affected / sqrt(n_affected),
    se_aquaculture_not_affected = sd_aquaculture_not_affected / sqrt(n_not_affected),
    
    # Compute standard error (SE) for salinity
    se_salinity_affected = sd_salinity_affected / sqrt(n_affected),
    se_salinity_not_affected = sd_salinity_not_affected / sqrt(n_not_affected),
    
    # Compute 95% Confidence Intervals for aquaculture
    upper_bound_aquaculture_affected = mean_aquaculture_affected + (1.96 * se_aquaculture_affected),
    lower_bound_aquaculture_affected = mean_aquaculture_affected - (1.96 * se_aquaculture_affected),
    
    upper_bound_aquaculture_not_affected = mean_aquaculture_not_affected + (1.96 * se_aquaculture_not_affected),
    lower_bound_aquaculture_not_affected = mean_aquaculture_not_affected - (1.96 * se_aquaculture_not_affected),
    
    # Compute 95% Confidence Intervals for salinity
    upper_bound_salinity_affected = mean_salinity_affected + (1.96 * se_salinity_affected),
    lower_bound_salinity_affected = mean_salinity_affected - (1.96 * se_salinity_affected),
    
    upper_bound_salinity_not_affected = mean_salinity_not_affected + (1.96 * se_salinity_not_affected),
    lower_bound_salinity_not_affected = mean_salinity_not_affected - (1.96 * se_salinity_not_affected)
  )

print(yearly_summary)


# Visualize aquaculture trends over time for affected and not affected villages
ggplot(yearly_summary) +
  # Add confidence bands using geom_ribbon (or geom_errorbar())
  geom_ribbon(aes(x = Year, ymin = lower_bound_aquaculture_affected, ymax = upper_bound_aquaculture_affected, fill = "Affected by Storms"), alpha = 0.2) +
  geom_ribbon(aes(x = Year, ymin = lower_bound_aquaculture_not_affected, ymax = upper_bound_aquaculture_not_affected, fill = "Not Affected by Storms"), alpha = 0.2) +
  
  # Lines for means
  geom_line(aes(x = Year, y = mean_aquaculture_affected, color = "Affected by Storms")) +
  geom_line(aes(x = Year, y = mean_aquaculture_not_affected, color = "Not Affected by Storms")) +
  
  # Points for means
  geom_point(aes(x = Year, y = mean_aquaculture_affected, color = "Affected by Storms")) +
  geom_point(aes(x = Year, y = mean_aquaculture_not_affected, color = "Not Affected by Storms")) +
  
  # Titles & labels
  labs(title = "Average Aquaculture Percentage by Year (Affected vs. Not Affected)",
       x = "Year", y = "Average Aquaculture (%)", color = "Storm Impact", fill = "Storm Impact") +
  scale_x_continuous(breaks = seq(min(yearly_summary$Year), max(yearly_summary$Year), by = 1)) +  # Ensure integer year labels
  
  # Custom colors for lines and fills
  scale_color_manual(values = c("Affected by Storms" = "blue", "Not Affected by Storms" = "orange")) +
  scale_fill_manual(values = c("Affected by Storms" = "blue", "Not Affected by Storms" = "orange")) +
  
  theme_minimal()

ggsave("outputs/Stormwise_aquaTrends_2013-25_allVillages.png", width = 10, height = 6, dpi = 300)

# Visualize salinity trends over time for affected and not affected villages
ggplot(yearly_summary, aes(x=factor(Year))) +
  geom_ribbon(aes(x = Year, ymin = lower_bound_salinity_affected, ymax = upper_bound_salinity_affected, fill = "Affected by Storms"), alpha = 0.2) +
  geom_ribbon(aes(x = Year, ymin = lower_bound_salinity_not_affected, ymax = upper_bound_salinity_not_affected, fill = "Not Affected by Storms"), alpha = 0.2) +
  geom_line(aes(x = Year, y = mean_salinity_affected, color = "Affected by Storms")) +
  geom_line(aes(x = Year, y = mean_salinity_not_affected, color = "Not Affected by Storms")) +
  geom_point(aes(x = Year, y = mean_salinity_affected, color = "Affected by Storms")) +
  geom_point(aes(x = Year, y = mean_salinity_not_affected, color = "Not Affected by Storms")) +
  scale_color_manual(values = c("Affected by Storms" = "blue", "Not Affected by Storms" = "orange")) +
  scale_fill_manual(values = c("Affected by Storms" = "blue", "Not Affected by Storms" = "orange")) +
  scale_x_continuous(breaks = seq(min(yearly_summary$Year), max(yearly_summary$Year), by = 1)) +  # Ensure integer year labels
  labs(title = "Average Saline Area Percentage by Year (Affected vs. Not Affected)",
       x = "Year", y = "Average Saline Area (%)", color = "Storm Impact", fill = "Storm Impact") +
  theme_minimal()
ggsave("outputs/Stormwise_salinityTrends_2013-25_allVillages.png", width = 10, height = 6, dpi = 300)




# NEAR VILLAGES (<60 km from the sea) 

# Summary stats for near and very near villages only 
summary_stats <- df_near %>%
  group_by(Year) %>%
  summarise(
    mean_aqua = mean(Aqua_perc, na.rm = TRUE),
    median_aqua = median(Aqua_perc, na.rm = TRUE),
    min_aqua = min(Aqua_perc, na.rm = TRUE),
    max_aqua = max(Aqua_perc, na.rm = TRUE),
    mean_saline = mean(Saline_perc, na.rm = TRUE),
    median_saline = median(Saline_perc, na.rm = TRUE),
    min_saline = min(Saline_perc, na.rm = TRUE),
    max_saline = max(Saline_perc, na.rm = TRUE)
  )

print(summary_stats)
# Finding: Overall trend shows increasing aquaculture and salinity in the near and very near villages over the decade



# near villages (<60km ) Divide for storm affected and unaffected
# Divide this for storm affected and unaffected
yearly_summary <- df_near %>%
  group_by(Year) %>%
  summarize(
    # Mean aquaculture % for affected and not affected villages
    mean_aquaculture_affected = mean(Aqua_perc[postSurge == 1], na.rm = TRUE),
    mean_aquaculture_not_affected = mean(Aqua_perc[postSurge == 0], na.rm = TRUE),
    
    # Standard deviation for aquaculture
    sd_aquaculture_affected = sd(Aqua_perc[postSurge == 1], na.rm = TRUE),
    sd_aquaculture_not_affected = sd(Aqua_perc[postSurge == 0], na.rm = TRUE),
    
    # Mean salinity % for affected and not affected villages
    mean_salinity_affected = mean(Saline_perc[postSurge == 1], na.rm = TRUE),
    mean_salinity_not_affected = mean(Saline_perc[postSurge == 0], na.rm = TRUE),
    
    # Standard deviation for salinity
    sd_salinity_affected = sd(Saline_perc[postSurge == 1], na.rm = TRUE),
    sd_salinity_not_affected = sd(Saline_perc[postSurge == 0], na.rm = TRUE),
    
    # Count of affected and not affected villages
    n_affected = sum(postSurge == 1, na.rm = TRUE),
    n_not_affected = sum(postSurge == 0, na.rm = TRUE),
    
    # Count of storm-affected villages
    storm_affected_villages = sum(postSurge, na.rm = TRUE),
    
    # Total number of unique villages
    n_villages = n_distinct(UniqueID)
  ) %>%
  mutate(
    # Compute standard error (SE) for aquaculture
    se_aquaculture_affected = sd_aquaculture_affected / sqrt(n_affected),
    se_aquaculture_not_affected = sd_aquaculture_not_affected / sqrt(n_not_affected),
    
    # Compute standard error (SE) for salinity
    se_salinity_affected = sd_salinity_affected / sqrt(n_affected),
    se_salinity_not_affected = sd_salinity_not_affected / sqrt(n_not_affected),
    
    # Compute 95% Confidence Intervals for aquaculture
    upper_bound_aquaculture_affected = mean_aquaculture_affected + (1.96 * se_aquaculture_affected),
    lower_bound_aquaculture_affected = mean_aquaculture_affected - (1.96 * se_aquaculture_affected),
    
    upper_bound_aquaculture_not_affected = mean_aquaculture_not_affected + (1.96 * se_aquaculture_not_affected),
    lower_bound_aquaculture_not_affected = mean_aquaculture_not_affected - (1.96 * se_aquaculture_not_affected),
    
    # Compute 95% Confidence Intervals for salinity
    upper_bound_salinity_affected = mean_salinity_affected + (1.96 * se_salinity_affected),
    lower_bound_salinity_affected = mean_salinity_affected - (1.96 * se_salinity_affected),
    
    upper_bound_salinity_not_affected = mean_salinity_not_affected + (1.96 * se_salinity_not_affected),
    lower_bound_salinity_not_affected = mean_salinity_not_affected - (1.96 * se_salinity_not_affected)
  )

print(yearly_summary)


# Visualize aquaculture trends over time for affected and not affected villages
ggplot(yearly_summary) +
  # Add confidence bands using geom_ribbon (or geom_errorbar())
  geom_ribbon(aes(x = Year, ymin = lower_bound_aquaculture_affected, ymax = upper_bound_aquaculture_affected, fill = "Affected by Storms"), alpha = 0.2) +
  geom_ribbon(aes(x = Year, ymin = lower_bound_aquaculture_not_affected, ymax = upper_bound_aquaculture_not_affected, fill = "Not Affected by Storms"), alpha = 0.2) +
  
  # Lines for means
  geom_line(aes(x = Year, y = mean_aquaculture_affected, color = "Affected by Storms")) +
  geom_line(aes(x = Year, y = mean_aquaculture_not_affected, color = "Not Affected by Storms")) +
  
  # Points for means
  geom_point(aes(x = Year, y = mean_aquaculture_affected, color = "Affected by Storms")) +
  geom_point(aes(x = Year, y = mean_aquaculture_not_affected, color = "Not Affected by Storms")) +
  
  # Titles & labels
  labs(title = "Average Aquaculture Percentage for Villages less than 60km from the sea",
       x = "Year", y = "Average Aquaculture (%)", color = "Storm Impact", fill = "Storm Impact") +
  scale_x_continuous(breaks = seq(min(yearly_summary$Year), max(yearly_summary$Year), by = 1)) +  # Ensure integer year labels
  
  # Custom colors for lines and fills
  scale_color_manual(values = c("Affected by Storms" = "blue", "Not Affected by Storms" = "orange")) +
  scale_fill_manual(values = c("Affected by Storms" = "blue", "Not Affected by Storms" = "orange")) +
  
  theme_minimal()
ggsave("outputs/Stormwise_aquaTrends_2013-25_nearVillages.png", width = 10, height = 6, dpi = 300)

# Visualize salinity trends over time for affected and not affected villages
ggplot(yearly_summary, aes(x=factor(Year))) +
  geom_ribbon(aes(x = Year, ymin = lower_bound_salinity_affected, ymax = upper_bound_salinity_affected, fill = "Affected by Storms"), alpha = 0.2) +
  geom_ribbon(aes(x = Year, ymin = lower_bound_salinity_not_affected, ymax = upper_bound_salinity_not_affected, fill = "Not Affected by Storms"), alpha = 0.2) +
  geom_line(aes(x = Year, y = mean_salinity_affected, color = "Affected by Storms")) +
  geom_line(aes(x = Year, y = mean_salinity_not_affected, color = "Not Affected by Storms")) +
  geom_point(aes(x = Year, y = mean_salinity_affected, color = "Affected by Storms")) +
  geom_point(aes(x = Year, y = mean_salinity_not_affected, color = "Not Affected by Storms")) +
  scale_color_manual(values = c("Affected by Storms" = "blue", "Not Affected by Storms" = "orange")) +
  scale_fill_manual(values = c("Affected by Storms" = "blue", "Not Affected by Storms" = "orange")) +
  scale_x_continuous(breaks = seq(min(yearly_summary$Year), max(yearly_summary$Year), by = 1)) +  # Ensure integer year labels
  labs(title = "Average Saline Area Percentage for Villages less than 60km from the sea",
       x = "Year", y = "Average Saline Area (%)", color = "Storm Impact", fill = "Storm Impact") +
  theme_minimal()
ggsave("outputs/Stormwise_salinityTrends_2013-25_nearVillages.png", width = 10, height = 6, dpi = 300)


# VERY NEAR VILLAGES (<30 km from the sea) 

# Summary stats for very near villages only 
summary_stats <- df_vnear %>%
  group_by(Year) %>%
  summarise(
    mean_aqua = mean(Aqua_perc, na.rm = TRUE),
    median_aqua = median(Aqua_perc, na.rm = TRUE),
    min_aqua = min(Aqua_perc, na.rm = TRUE),
    max_aqua = max(Aqua_perc, na.rm = TRUE),
    mean_saline = mean(Saline_perc, na.rm = TRUE),
    median_saline = median(Saline_perc, na.rm = TRUE),
    min_saline = min(Saline_perc, na.rm = TRUE),
    max_saline = max(Saline_perc, na.rm = TRUE)
  )

print(summary_stats)
# Finding: Overall trend shows increasing aquaculture and salinity in the near and very near villages over the decade



# very near villages (<30km ) Divide for storm affected and unaffected
yearly_summary <- df_vnear %>%
  group_by(Year) %>%
  summarize(
    # Mean aquaculture % for affected and not affected villages
    mean_aquaculture_affected = mean(Aqua_perc[postSurge == 1], na.rm = TRUE),
    mean_aquaculture_not_affected = mean(Aqua_perc[postSurge == 0], na.rm = TRUE),
    
    # Standard deviation for aquaculture
    sd_aquaculture_affected = sd(Aqua_perc[postSurge == 1], na.rm = TRUE),
    sd_aquaculture_not_affected = sd(Aqua_perc[postSurge == 0], na.rm = TRUE),
    
    # Mean salinity % for affected and not affected villages
    mean_salinity_affected = mean(Saline_perc[postSurge == 1], na.rm = TRUE),
    mean_salinity_not_affected = mean(Saline_perc[postSurge == 0], na.rm = TRUE),
    
    # Standard deviation for salinity
    sd_salinity_affected = sd(Saline_perc[postSurge == 1], na.rm = TRUE),
    sd_salinity_not_affected = sd(Saline_perc[postSurge == 0], na.rm = TRUE),
    
    # Count of affected and not affected villages
    n_affected = sum(postSurge == 1, na.rm = TRUE),
    n_not_affected = sum(postSurge == 0, na.rm = TRUE),
    
    # Count of storm-affected villages
    storm_affected_villages = sum(postSurge, na.rm = TRUE),
    
    # Total number of unique villages
    n_villages = n_distinct(UniqueID)
  ) %>%
  mutate(
    # Compute standard error (SE) for aquaculture
    se_aquaculture_affected = sd_aquaculture_affected / sqrt(n_affected),
    se_aquaculture_not_affected = sd_aquaculture_not_affected / sqrt(n_not_affected),
    
    # Compute standard error (SE) for salinity
    se_salinity_affected = sd_salinity_affected / sqrt(n_affected),
    se_salinity_not_affected = sd_salinity_not_affected / sqrt(n_not_affected),
    
    # Compute 95% Confidence Intervals for aquaculture
    upper_bound_aquaculture_affected = mean_aquaculture_affected + (1.96 * se_aquaculture_affected),
    lower_bound_aquaculture_affected = mean_aquaculture_affected - (1.96 * se_aquaculture_affected),
    
    upper_bound_aquaculture_not_affected = mean_aquaculture_not_affected + (1.96 * se_aquaculture_not_affected),
    lower_bound_aquaculture_not_affected = mean_aquaculture_not_affected - (1.96 * se_aquaculture_not_affected),
    
    # Compute 95% Confidence Intervals for salinity
    upper_bound_salinity_affected = mean_salinity_affected + (1.96 * se_salinity_affected),
    lower_bound_salinity_affected = mean_salinity_affected - (1.96 * se_salinity_affected),
    
    upper_bound_salinity_not_affected = mean_salinity_not_affected + (1.96 * se_salinity_not_affected),
    lower_bound_salinity_not_affected = mean_salinity_not_affected - (1.96 * se_salinity_not_affected)
  )

print(yearly_summary)


# Visualize aquaculture trends over time for affected and not affected villages
ggplot(yearly_summary) +
  # Add confidence bands using geom_ribbon (or geom_errorbar())
  geom_ribbon(aes(x = Year, ymin = lower_bound_aquaculture_affected, ymax = upper_bound_aquaculture_affected, fill = "Affected by Storms"), alpha = 0.2) +
  geom_ribbon(aes(x = Year, ymin = lower_bound_aquaculture_not_affected, ymax = upper_bound_aquaculture_not_affected, fill = "Not Affected by Storms"), alpha = 0.2) +
  
  # Lines for means
  geom_line(aes(x = Year, y = mean_aquaculture_affected, color = "Affected by Storms")) +
  geom_line(aes(x = Year, y = mean_aquaculture_not_affected, color = "Not Affected by Storms")) +
  
  # Points for means
  geom_point(aes(x = Year, y = mean_aquaculture_affected, color = "Affected by Storms")) +
  geom_point(aes(x = Year, y = mean_aquaculture_not_affected, color = "Not Affected by Storms")) +
  
  # Titles & labels
  labs(title = "Average Aquaculture Percentage for Villages less than 30km from the sea",
       x = "Year", y = "Average Aquaculture (%)", color = "Storm Impact", fill = "Storm Impact") +
  scale_x_continuous(breaks = seq(min(yearly_summary$Year), max(yearly_summary$Year), by = 1)) +  # Ensure integer year labels
  
  # Custom colors for lines and fills
  scale_color_manual(values = c("Affected by Storms" = "blue", "Not Affected by Storms" = "orange")) +
  scale_fill_manual(values = c("Affected by Storms" = "blue", "Not Affected by Storms" = "orange")) +
  
  theme_minimal()
ggsave("outputs/Stormwise_aquaTrends_2013-25_verynearVillages.png", width = 10, height = 6, dpi = 300)


# Visualize salinity trends over time for affected and not affected villages
ggplot(yearly_summary, aes(x=factor(Year))) +
  geom_ribbon(aes(x = Year, ymin = lower_bound_salinity_affected, ymax = upper_bound_salinity_affected, fill = "Affected by Storms"), alpha = 0.2) +
  geom_ribbon(aes(x = Year, ymin = lower_bound_salinity_not_affected, ymax = upper_bound_salinity_not_affected, fill = "Not Affected by Storms"), alpha = 0.2) +
  geom_line(aes(x = Year, y = mean_salinity_affected, color = "Affected by Storms")) +
  geom_line(aes(x = Year, y = mean_salinity_not_affected, color = "Not Affected by Storms")) +
  geom_point(aes(x = Year, y = mean_salinity_affected, color = "Affected by Storms")) +
  geom_point(aes(x = Year, y = mean_salinity_not_affected, color = "Not Affected by Storms")) +
  scale_color_manual(values = c("Affected by Storms" = "blue", "Not Affected by Storms" = "orange")) +
  scale_fill_manual(values = c("Affected by Storms" = "blue", "Not Affected by Storms" = "orange")) +
  scale_x_continuous(breaks = seq(min(yearly_summary$Year), max(yearly_summary$Year), by = 1)) +  # Ensure integer year labels
  labs(title = "Average Saline Area Percentage for Villages less than 30km from the sea",
       x = "Year", y = "Average Saline Area (%)", color = "Storm Impact", fill = "Storm Impact") +
  theme_minimal()

ggsave("outputs/Stormwise_salinityTrends_2013-25_verynearVillages.png", width = 10, height = 6, dpi = 300)



# 3. Time series for Mapping all villages
storm_effect <- df %>%
  group_by(postSurge) %>%  # Binary: 0 = No, 1 = Yes
  summarise(
    mean_aquaculture = mean(Aqua_perc, na.rm = TRUE),
    mean_saline = mean(Saline_perc, na.rm = TRUE),
    mean_dry_aquaculture = mean(DryAqua_perc, na.rm = TRUE)
  )

print(storm_effect)
# Finding: Storm affected areas have higher aquaculture, dry aquaculture as well as salinity prevelance 


# Compute mean aquaculture percentage per year for storm-affected and unaffected villages
df_mean <- df %>%
  group_by(Year, postSurge) %>%
  summarise(mean_aqua = mean(Aqua_perc, na.rm = TRUE), .groups = "drop")

# Plot with individual village trends and mean trend lines
ggplot(df, aes(x = Year, y = Aqua_perc, group = UniqueID, color = factor(postSurge))) +
  geom_line(alpha = 0.3) +  # Individual village trends
  geom_line(data = df_mean, aes(x = Year, y = mean_aqua, group = postSurge, color = factor(postSurge)), 
            size = 1.2) +  # Mean trend line
  labs(x = "Year", y = "% Aquaculture", 
       title = "Aquaculture Trends Over Time by Storm Affected and Unaffected Areas", 
       color = "Storm Surge Impact") +
  coord_cartesian(ylim = c(0, 5)) +  # Zoom into y-axis range (adjust)
  theme_minimal()



# For near and very near villages only 
storm_effect <- df_near %>%
  group_by(postSurge) %>%  # Binary: 0 = No, 1 = Yes
  summarise(
    mean_aquaculture = mean(Aqua_perc, na.rm = TRUE),
    mean_saline = mean(Saline_perc, na.rm = TRUE),
    mean_dry_aquaculture = mean(DryAqua_perc, na.rm = TRUE)
  )

print(storm_effect)
# Finding: Storm affected areas have higher aquaculture, dry aquaculture as well as salinity prevelance as compared to unaffected areas, despite being similarly close to the sea. 

# Time series for villages near and very near to the sea
# Compute for df_near villages - mean aquaculture percentage per year for storm-affected and unaffected villages 
df_mean <- df_near %>%
  group_by(Year, postSurge) %>%
  summarise(mean_aqua = mean(Aqua_perc, na.rm = TRUE), .groups = "drop")

# Plot with individual village trends and mean trend lines
ggplot(df_near, aes(x = Year, y = Aqua_perc, group = UniqueID, color = factor(postSurge))) +
  geom_line(alpha = 0.3) +  # Individual village trends
  geom_line(data = df_mean, aes(x = Year, y = mean_aqua, group = postSurge, color = factor(postSurge)), 
            size = 1.2) +  # Mean trend line
  labs(x = "Year", y = "% Aquaculture", 
       title = "Aquaculture Trends Over Time by Storm Affected and Unaffected Areas close to the sea (within 60km)", 
       color = "Storm Surge Impact") +
  coord_cartesian(ylim = c(0, 5)) +  # Zoom into y-axis range (adjust)
  theme_minimal()



# For very near villages only
storm_effect <- df_vnear %>%
  group_by(postSurge) %>%  # Binary: 0 = No, 1 = Yes
  summarise(
    mean_aquaculture = mean(Aqua_perc, na.rm = TRUE),
    mean_saline = mean(Saline_perc, na.rm = TRUE),
    mean_dry_aquaculture = mean(DryAqua_perc, na.rm = TRUE)
  )

print(storm_effect)
# Finding: Storm affected areas have higher aquaculture, dry aquaculture as well as salinity prevelance as compared to unaffected areas, despite being similarly close (<30 km) to the sea.

# Compute for df_vnear villages - mean aquaculture percentage per year for storm-affected and unaffected villages
df_mean <- df_vnear %>%
  group_by(Year, postSurge) %>%
  summarise(mean_aqua = mean(Aqua_perc, na.rm = TRUE), .groups = "drop")

# Plot with individual village trends and mean trend lines
ggplot(df_vnear, aes(x = Year, y = Aqua_perc, group = UniqueID, color = factor(postSurge))) +
  geom_line(alpha = 0.3) +  # Individual village trends
  geom_line(data = df_mean, aes(x = Year, y = mean_aqua, group = postSurge, color = factor(postSurge)),
            size = 1.2) +  # Mean trend line
  labs(x = "Year", y = "% Aquaculture",
       title = "Aquaculture Trends Over Time by Storm Affected and Unaffected Areas close to the sea (within 30km)",
       color = "Storm Surge Impact") +
  coord_cartesian(ylim = c(0, 5)) +  # Zoom into y-axis range (adjust)
  theme_minimal()





# 4. Aquaculture by both saline and storm affected areas 
# Step 1: Use binary variable for saline places - above and below the median (too many extreme values)
summary(df$Saline)
table(df$Saline_Storm_Category)

# Step 2: Summarize aquaculture trends
yearly_summary <- df %>%
  group_by(Year, Saline_Storm_Category) %>%
  summarize(
    mean_aquaculture = mean(Aqua_perc, na.rm = TRUE),
    sd_aquaculture = sd(Aqua_perc, na.rm = TRUE),
    upper_bound = mean_aquaculture + sd_aquaculture,
    lower_bound = mean_aquaculture - sd_aquaculture,
    .groups = "drop"
  )

# Step 3: Plot aquaculture trends over time

# Define custom colors
category_colors <- c(
  "Flooded & High Saline" = "blue",
  "Non-Flooded & High Saline" = "orange",
  "Flooded & Low Saline" = "darkgreen",
  "Non-Flooded & Low Saline" = "grey"
)

ggplot(yearly_summary, aes(x = Year, y = mean_aquaculture, color = Saline_Storm_Category, fill = Saline_Storm_Category)) +
  #geom_ribbon(aes(ymin = lower_bound, ymax = upper_bound), alpha = 0.2) +  # Error bands
  geom_line(size = 1) +
  geom_point(size = 2) +
  labs(title = "Aquaculture Trends Over Time by Salinity & Storm Impact",
       x = "Year", y = "Average Aquaculture (%)", color = "Category", fill = "Category") +
  scale_x_continuous(breaks = seq(min(yearly_summary$Year), max(yearly_summary$Year), by = 1)) +  # Ensure integer year labels
  scale_color_manual(values = category_colors) +  # Apply custom colors
  scale_fill_manual(values = category_colors) +  # Match fill colors to lines
  theme_minimal()

# Save plot
ggsave("outputs/Aquaculture_Trends_Salinity_Storm.png", width = 10, height = 6, dpi = 300)



# 5. Correlations Matrix
# All village
cor_matrix <- df %>%
  select(Aqua_perc, Saline_perc, DryAqua_perc, postSurge) %>%
  cor(use = "complete.obs")

corrplot(cor_matrix, method = "color", type = "lower", tl.cex = 0.8)
# Observation: Dry Aquaculture is highly correlated with salinity 

# Near and very near villages
cor_matrix <- df_near %>%
  select(Aqua_perc, Saline_perc, DryAqua_perc, postSurge) %>%
  cor(use = "complete.obs")

corrplot(cor_matrix, method = "color", type = "lower", tl.cex = 0.8)
# Observation: Dry Aquaculture is even more correlated with salinity 


# Very near villages
cor_matrix <- df_vnear %>%
  select(Aqua_perc, Saline_perc, DryAqua_perc, postSurge) %>%
  cor(use = "complete.obs")

corrplot(cor_matrix, method = "color", type = "lower", tl.cex = 0.8)
# Observation: Similar to near and very near villages



# 6. Lagged Correlation (Salinity t-1 vs. Aquaculture t)
df_lag <- df %>%
  group_by(UniqueID) %>%
  mutate(Lag_Saline = lag(Saline_perc))  # Lag salinity by one year

cor.test(df_lag$Lag_Saline, df_lag$Aqua_perc, use = "complete.obs")
# Interpretation: While the correlation is statistically significant (p < 0.05), the effect size is quite small (r ≈ 0.16). This means:
# There is a real relationship between the variables (not due to chance)
# But the relationship explains only about 2.6% of the variance (r² ≈ 0.026)
# The large sample size (n > 368,000) is likely why even this small correlation is highly significant

ggplot(df_lag, aes(x = Lag_Saline, y = Aqua_perc)) +
  geom_point(alpha = 0.5) +
  geom_smooth(method = "lm", se = FALSE, color = "red") +
  labs(x = "Salinity in Previous Year", y = "Aquaculture in Current Year",
       title = "Impact of Previous Year's Salinity on Aquaculture") +
  theme_minimal()



# 7. Lagged Correlation (Aquaculture t-1 vs. Salinity t)
df_lag <- df %>%
  group_by(UniqueID) %>%
  mutate(Lag_Aqua = lag(Aqua_perc))  # Lag salinity by one year

cor.test(df_lag$Lag_Aqua, df_lag$Saline_perc, use = "complete.obs")
# Interpretation: this result is statistically significant (p < 0.05) but shows a small effect size (r ≈ 0.15). This means:
# There is a real relationship between lag in aqua and saline percentage
# The relationship explains only about 2.2% of the variance (r² ≈ 0.022)
# The large sample size (n > 368,000) is likely why this small correlation is highly significant

ggplot(df_lag, aes(x = Lag_Aqua, y = Saline_perc)) +
  geom_point(alpha = 0.5) +
  geom_smooth(method = "lm", se = FALSE, color = "red") +
  labs(x = "Aquaculture in Previous Year", y = "Salinity in Current Year",
       title = "Impact of Previous Year's Aquaculture on Salinity") +
  theme_minimal()
# Observation: The linear correlation indicates a positive relationship between aquaculture in the previous year and salinity in the current year


# 8. Persistence of salinity over 5 years and relationship with aquaculture
# Test not just lag by one year, but average salinity in the last 5 years with aquaculture 
head(df)

df_lagged <- df %>%
  group_by(UniqueID) %>%
  arrange(Year) %>%
  mutate(
    Salinity_t1 = lag(Saline_perc, 1),
    Salinity_t2 = lag(Saline_perc, 2),
    Salinity_t3 = lag(Saline_perc, 3),
    Salinity_t4 = lag(Saline_perc, 4),
    Salinity_t5 = lag(Saline_perc, 5),
    Avg_SalinityPerc_Last5 = rowMeans(cbind(Salinity_t1, Salinity_t2, Salinity_t3, Salinity_t4, Salinity_t5), na.rm = TRUE)
  )

cor.test(df_lagged$Avg_SalinityPerc_Last5, df_lagged$Aqua_perc, use = "complete.obs")
# Interpretation: The relationship is statistically significant, but significance does not imply strength—the correlation is still weak (0.1957). 
# This suggests that higher past salinity might be slightly linked to increased aquaculture, but salinity alone is not a strong predictor.
# Other factors (e.g., geography, storms, policy, land use, economic incentives) might also be influencing aquaculture expansion.


# 9. Check for non-linearity 
ggplot(df_lagged, aes(x = Avg_SalinityPerc_Last5, y = Aqua_perc)) +
  geom_point(alpha = 0.2) +
  geom_smooth(method = "loess", color = "red") +
  theme_minimal() +
  labs(title = "Relationship between Past Salinity (5-Year Avg) and Aquaculture")


# Next step is to run multiple regression models to account for other variables and isolate the effects of salinity and storm surge 
