# Compile distance from the sea 

library(dplyr)

getwd()

file.choose()
OD_dist <- read.csv("outputs/OD_Village_Distance_To_Sea.csv")
AP_dist <- read.csv("outputs/AP_Village_Distance_To_Sea.csv")
TN_dist <- read.csv("outputs/TN_Village_Distance_To_Sea.csv")

head(OD_dist)
head(AP_dist)
head(TN_dist)

ALL_dist <- bind_rows(OD_dist, AP_dist, TN_dist)

summary(ALL_dist) # should be 31088
table(duplicated(ALL_dist$UniqueID))  # TRUE means duplicates exist
# select the first value only
ALL_dist <- ALL_dist %>% distinct(UniqueID, .keep_all = TRUE)

summary(ALL_dist)

write.csv(ALL_dist, "outputs/all_SeaDistances.csv")
