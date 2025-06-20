# Compare wealth data for the select districts 
# DHS Wealth Data 

library(ipumsr)
library(readr)
library(dplyr)
library(ggplot2)


wealth_layout <- fwf_widths(
  widths = c(5, 5, 3, 4, 21, 14, 11, 10, 15, 2, 9, 10, 8, 2, 5, 1, 1, 12, 1, 8),
  col_names = c(
    "SAMPLE", "SAMPLESTR", "COUNTRY", "YEAR", "IDHSHID", "DHSID", "IDHSPSU", "IDHSSTRATA",
    "HHID", "HHLINENO", "HHWEIGHT", "POPWT_HH", "HHMWEIGHT", "GEO_IA1992_2019", 
    "GEOALT_IA2015_2019", "HHRESIDENT", "WEALTHQHH", "WEALTHSHH", 
    "WEALTHQRURBHH", "WEALTHSRURBHH"
  )
)

wealth <- read_fwf("C:\\Users\\garim\\Documents\\GitHub_LocalRepository\\aquaculture_patterns\\data\\DHS India\\Wealth\\idhs_00002.dat\\idhs_00002.dat", wealth_layout)

head(wealth)

# GEO_IA1992_2019 == 21 for OD, 28 for AP, 33 for TN 
# GEOALT_IA2015_2019 for district codes (note - only available for 2015 and 2019, rest NA) == 21008, 21009, 21010, 21011,21018,21019,28001,28002,28003,28004,28005,28006,28007,28008,28009,33001,33003,33006,33016,33017,33018,33019,33020,33025,33026,33027,33028
# WEALTHQHH (quintile: 1 = poorest, 5 = richest)

print(unique(wealth$YEAR))

# Filter dataset to include only necessary entries
district_codes <- c(21008, 21009, 21010, 21011, 21018, 21019,
                    28001, 28002, 28003, 28004, 28005, 28006, 28007, 28008, 28009,
                    33001, 33003, 33006, 33016, 33017, 33018, 33019, 33020,
                    33025, 33026, 33027, 33028)

state_codes <- c(21, 28, 33)
years <- c(1992, 1998, 2005, 2015, 2019)


wealth_filtered <- wealth %>%
  filter(YEAR %in% years,
         GEO_IA1992_2019 %in% state_codes) %>%
  # For 2015/2019 restrict to selected districts
  filter(
    (YEAR %in% c(2015, 2019) & GEOALT_IA2015_2019 %in% district_codes) |
      (YEAR %in% c(1992, 1998, 2005))  # earlier years: accept all districts (if not available)
  )


# Add labels for clarity
wealth_filtered <- wealth_filtered %>%
  mutate(
    State = case_when(
      GEO_IA1992_2019 == 21 ~ "OD",
      GEO_IA1992_2019 == 28 ~ "AP",
      GEO_IA1992_2019 == 33 ~ "TN"
    )
  )

table(wealth$WEALTHQHH, wealth$YEAR)


# State wealth averages
state_avg_wealth <- wealth_filtered %>%
  group_by(YEAR, State) %>%
  summarise(
    n = sum(!is.na(WEALTHQHH)),
    mean_quintile = mean(as.numeric(WEALTHQHH), na.rm = TRUE),
    sd_quintile = sd(as.numeric(WEALTHQHH), na.rm = TRUE),
    se_quintile = sd_quintile / sqrt(n),
    ci_lower = mean_quintile - 1.96 * se_quintile,
    ci_upper = mean_quintile + 1.96 * se_quintile,
    .groups = "drop"
  )

print(state_avg_wealth)

ggplot(state_avg_wealth, aes(x = YEAR, y = mean_quintile, color = State)) +
  geom_line(size = 1.2, aes(group = State)) +
  geom_point(size = 2) +
  geom_ribbon(aes(ymin = ci_lower, ymax = ci_upper, fill = State, group = State), alpha = 0.2, color = NA) +
  scale_y_continuous(breaks = 1:5, limits = c(1, 5)) +
  scale_x_continuous(breaks = c(1992, 1998, 2005, 2015, 2019)) +
  labs(
    title = "Average Wealth Quintile Over Time by State (with 95% CI)",
    x = "Survey Year",
    y = "Average Wealth Quintile",
    color = "State",
    fill = "State"
  ) +
  theme_minimal(base_size = 14)+
  theme(
    legend.position = "bottom",
    legend.direction = "horizontal",
    legend.box = "horizontal"
  )


# District wealth averages (2015 and 2019 only) 

wealth_districts <- wealth_filtered %>%
  filter(
    YEAR %in% c(2015, 2019),
    as.numeric(GEOALT_IA2015_2019) %in% district_codes
  ) %>%
  mutate(
    District = as.factor(GEOALT_IA2015_2019),
    State = case_when(
      as.numeric(GEOALT_IA2015_2019) %in% c(21008, 21009, 21010, 21011, 21018, 21019) ~ "OD",
      as.numeric(GEOALT_IA2015_2019) %in% c(28001, 28002, 28003, 28004, 28005, 28006, 28007, 28008, 28009) ~ "AP",
      as.numeric(GEOALT_IA2015_2019) %in% c(33001, 33003, 33006, 33016, 33017, 33018, 33019, 33020, 33025, 33026, 33027, 33028) ~ "TN",
      TRUE ~ NA_character_
    )
  )

district_avg_wealth_ci <- wealth_districts %>%
  select(IDHSHID, District,State, YEAR, WEALTHQHH) %>%
  group_by(YEAR, District, State) %>%
  summarise(
    n = sum(!is.na(WEALTHQHH)),
    mean_quintile = mean(as.numeric(WEALTHQHH), na.rm = TRUE),
    sd_quintile = sd(as.numeric(WEALTHQHH), na.rm = TRUE),
    se_quintile = sd_quintile / sqrt(n),
    ci_lower = mean_quintile - 1.96 * se_quintile,
    ci_upper = mean_quintile + 1.96 * se_quintile,
    .groups = "drop"
  )

ggplot(district_avg_wealth_ci, aes(x = YEAR, y = mean_quintile, group = District)) +
  geom_line(aes(color = State), size = 1) +
  geom_point(aes(color = State), size = 2) +
  geom_ribbon(aes(ymin = ci_lower, ymax = ci_upper, fill = State), alpha = 0.2, color = NA) +
  scale_color_manual(values = c("OD" = "darkgreen", "AP" = "orange", "TN" = "blue")) +
  scale_fill_manual(values = c("OD" = "darkgreen", "AP" = "orange", "TN" = "blue")) +
  scale_y_continuous(breaks = 1:5, limits = c(1, 5)) +
  scale_x_continuous(breaks = c(2015, 2019)) +
  labs(
    title = "Average Wealth Quintile by District and State (2015 & 2019)",
    x = "Year",
    y = "Mean Wealth Quintile",
    color = "State",
    fill = "State"
  ) +
  theme_minimal(base_size = 13)
