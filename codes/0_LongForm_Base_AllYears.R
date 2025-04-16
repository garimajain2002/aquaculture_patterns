# Make a long form of all villages with all years as the base to bind other data to

unique_ids <- read.csv("data/All_uniqueIDs.csv")

unique_ids <- unique_ids %>%
  select(UniqueID, state_code, district_code) %>%
  distinct()

head(unique_ids)

summary(unique_ids)


# Create a dataframe for the years 1990-2025
years <- 1990:2025


# Generate all combinations of UniqueID and historical years
villageYears <- unique_ids %>%
  crossing(year = years) %>%
  # Keep only state_code and district_code columns from unique_ids
  select(UniqueID, year, state_code, district_code)

head(villageYears)
summary(villageYears)

write.csv(villageYears, "data/All_villageYears.csv")
