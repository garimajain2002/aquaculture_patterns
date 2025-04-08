#---- Calculate change in aqua and saline areas by year --------------
library(ggplot2)
library(tidyr)

setwd("C:/Users/garim/Documents/GitHub_LocalRepository/aquaculture_patterns")
getwd()

# ========== JAGATSINGHPUR TREND 

df_1995 <- read.csv(unz("outputs/1995_JSP_predicted_ECAqua_df.zip","1995_JSP_predicted_ECAqua_df.csv" ))
df_2001 <- read.csv(unz("outputs/2001_JSP_predicted_ECAqua_df.zip","2001_JSP_predicted_ECAqua_df.csv" ))
df_2005 <- read.csv(unz("outputs/2005_JSP_predicted_ECAqua_df.zip","2005_JSP_predicted_ECAqua_df.csv" ))
df_2010 <- read.csv(unz("outputs/2010_JSP_predicted_ECAqua_df.zip","2010_JSP_predicted_ECAqua_df.csv" ))
df_2014 <- read.csv(unz("outputs/2014_JSP_predicted_ECAqua_df.zip","2014_JSP_predicted_ECAqua_df.csv" ))
df_2017 <- read.csv(unz("outputs/2017_JSP_predicted_ECAqua_df.zip","2017_JSP_predicted_ECAqua_df.csv" ))
df_2021 <- read.csv(unz("outputs/2021_JSP_predicted_ECAqua_df.zip","2021_JSP_predicted_ECAqua_df.csv" ))
df_2024 <- read.csv(unz("outputs/2024_JSP_predicted_ECAqua_df.zip","2024_JSP_predicted_ECAqua_df.csv" ))

# Calculate area under aquaculture and salinity - convert from 30mx30m cell to sq. km 
aqua_1995 <- sum(df_1995$predicted_EC==2, na.rm=TRUE)*0.0009
saline_1995 <- sum(df_1995$predicted_EC==1, na.rm=TRUE)*0.0009

aqua_2001 <- sum(df_2001$predicted_EC==2, na.rm=TRUE)*0.0009
saline_2001 <- sum(df_2001$predicted_EC==1, na.rm=TRUE)*0.0009

aqua_2005 <- sum(df_2005$predicted_EC==2, na.rm=TRUE)*0.0009
saline_2005 <- sum(df_2005$predicted_EC==1, na.rm=TRUE)*0.0009

aqua_2010 <- sum(df_2010$predicted_EC==2, na.rm=TRUE)*0.0009
saline_2010 <- sum(df_2010$predicted_EC==1, na.rm=TRUE)*0.0009

aqua_2014 <- sum(df_2014$predicted_EC==2, na.rm=TRUE)*0.0009
saline_2014 <- sum(df_2014$predicted_EC==1, na.rm=TRUE)*0.0009

aqua_2017 <- sum(df_2017$predicted_EC==2, na.rm=TRUE)*0.0009
saline_2017 <- sum(df_2017$predicted_EC==1, na.rm=TRUE)*0.0009

aqua_2021 <- sum(df_2021$predicted_EC==2, na.rm=TRUE)*0.0009
saline_2021 <- sum(df_2021$predicted_EC==1, na.rm=TRUE)*0.0009

aqua_2024 <- sum(df_2024$predicted_EC==2, na.rm=TRUE)*0.0009
saline_2024 <- sum(df_2024$predicted_EC==1, na.rm=TRUE)*0.0009


combined_df <- data.frame(
  Year = c(1995, 2001, 2005, 2010, 2014, 2017, 2021, 2024),
  Aquaculture = c(aqua_1995, aqua_2001, aqua_2005, aqua_2010, aqua_2014, aqua_2017, aqua_2021, aqua_2024),
  Saline_Soil = c(saline_1995, saline_2001, saline_2005, saline_2010, saline_2014, saline_2017, saline_2021, saline_2024)
)

print(combined_df)

# Reshape the data from wide to long format
long_df <- pivot_longer(combined_df, 
                        cols = c(Aquaculture, Saline_Soil),
                        names_to = "Category",
                        values_to = "Count")

# Create the plot
ggplot(long_df, aes(x = Year, y = Count, color = Category)) +
  geom_line(linewidth = 1) +
  geom_smooth(method = "lm", linetype = "dotted", se = FALSE) +  
  scale_color_manual(values = c("Aquaculture" = "blue", "Saline_Soil" = "orange")) +
  theme_minimal() +
  labs(title = "Aquaculture and Saline Soil (1995-2024): Jagatsinghpur District",
       y = "Area in Sq.km.",
       x = "Year") +
  theme(legend.title = element_blank())





# ========== STUDY AREA TREND 

df_1995 <- read.csv("outputs/1995_SA_predicted_ECAqua_df.csv" )
df_2001 <- read.csv("outputs/2001_SA_predicted_ECAqua_df.csv" )
df_2005 <- read.csv("outputs/2005_SA_predicted_ECAqua_df.csv" )
df_2010 <- read.csv("outputs/2010_SA_predicted_ECAqua_df.csv" )
df_2014 <- read.csv("outputs/2014_SA_predicted_ECAqua_df.csv" )
df_2017 <- read.csv("outputs/2017_SA_predicted_ECAqua_df.csv" )
df_2021 <- read.csv("outputs/2021_SA_predicted_ECAqua_df.csv" )
df_2024 <- read.csv("outputs/2024_SA_predicted_ECAqua_df.csv" )

# Calculate area under quaculture and salinity - convert from 30mx30m cell to sq. km 
aqua_1995 <- sum(df_1995$predicted_EC==2, na.rm=TRUE)*0.0009
saline_1995 <- sum(df_1995$predicted_EC==1, na.rm=TRUE)*0.0009

aqua_2001 <- sum(df_2001$predicted_EC==2, na.rm=TRUE)*0.0009
saline_2001 <- sum(df_2001$predicted_EC==1, na.rm=TRUE)*0.0009

aqua_2005 <- sum(df_2005$predicted_EC==2, na.rm=TRUE)*0.0009
saline_2005 <- sum(df_2005$predicted_EC==1, na.rm=TRUE)*0.0009

aqua_2010 <- sum(df_2010$predicted_EC==2, na.rm=TRUE)*0.0009
saline_2010 <- sum(df_2010$predicted_EC==1, na.rm=TRUE)*0.0009

aqua_2014 <- sum(df_2014$predicted_EC==2, na.rm=TRUE)*0.0009
saline_2014 <- sum(df_2014$predicted_EC==1, na.rm=TRUE)*0.0009

aqua_2017 <- sum(df_2017$predicted_EC==2, na.rm=TRUE)*0.0009
saline_2017 <- sum(df_2017$predicted_EC==1, na.rm=TRUE)*0.0009

aqua_2021 <- sum(df_2021$predicted_EC==2, na.rm=TRUE)*0.0009
saline_2021 <- sum(df_2021$predicted_EC==1, na.rm=TRUE)*0.0009

aqua_2024 <- sum(df_2024$predicted_EC==2, na.rm=TRUE)*0.0009
saline_2024 <- sum(df_2024$predicted_EC==1, na.rm=TRUE)*0.0009


combined_df <- data.frame(
  Year = c(1995, 2001, 2005, 2010, 2014, 2017, 2021, 2024),
  Aquaculture = c(aqua_1995, aqua_2001, aqua_2005, aqua_2010, aqua_2014, aqua_2017, aqua_2021, aqua_2024),
  Saline_Soil = c(saline_1995, saline_2001, saline_2005, saline_2010, saline_2014, saline_2017, saline_2021, saline_2024)
)

print(combined_df)

# Reshape the data from wide to long format
long_df <- pivot_longer(combined_df, 
                        cols = c(Aquaculture, Saline_Soil),
                        names_to = "Category",
                        values_to = "Count")

# Create the plot
ggplot(long_df, aes(x = Year, y = Count, color = Category)) +
  geom_line(linewidth = 1) +
  geom_smooth(method = "lm", linetype = "dotted", se = FALSE) +  
  scale_color_manual(values = c("Aquaculture" = "blue", "Saline_Soil" = "orange")) +
  theme_minimal() +
  labs(title = "Aquaculture and Saline Soil (1995-2024): Study Area",
       y = "Area in Sq.km.",
       x = "Year") +
  theme(legend.title = element_blank())
