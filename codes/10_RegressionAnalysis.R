# Aquaculture Analysis: Panel Data Regression Models
# This script analyzes the relationship between aquaculture development,
# storm surges, salinity, and other environmental factors

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

data <- read.csv("outputs/aqua_salinity_surge_2013-2025.csv")

head(data)
summary(data)

# Check for missing values
missing_values <- colSums(is.na(data))
print(missing_values)
# Note: >150k are missing pct change. 
# ! Check why and fix, or don't use that variable. 

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
    avg_salinity_5yr = rollapply(Saline_perc, 
                                 width = 5, 
                                 FUN = mean, 
                                 fill = NA, 
                                 align = "right", 
                                 na.rm = TRUE, 
                                 partial = TRUE)  # Allow smaller windows when data is missing
  ) %>%
  ungroup()

summary(data$avg_salinity_5yr)


# Create interaction term for ease
data$salinity_storm_interaction <- data$Saline_perc * data$postSurge
summary(data)

# Some villages have multiple boundaries but same UniqueID and Year combination. Take their average value 
duplicate_check <- data %>%
  group_by(UniqueID, Year) %>%
  summarise(n = n(), .groups = "drop") %>%
  filter(n > 1)

print(duplicate_check)

data_agg <- data %>%
  group_by(UniqueID, Year) %>%
  summarize(
    Aqua_ha = mean(Aqua_ha, na.rm = TRUE),
    Aqua_perc = mean(Aqua_perc, na.rm = TRUE),
    Saline_perc = mean(Saline_perc, na.rm = TRUE),
    postSurge = max(postSurge, na.rm = TRUE), # A village is affected if any record says so
    salinity_storm_interaction <- Saline_perc * postSurge,
    DEM_avg = mean(DEM_avg, na.rm = TRUE), 
    NEAR_DIST = mean(NEAR_DIST, na.rm = TRUE),
    
    .groups = "drop"
  )

summary(data_agg)

# Set up panel data structure
pdata <- pdata.frame(data_agg, index = c("UniqueID", "Year"))



# ----------------------
# REGRESSION MODELS
# ----------------------
# RQ: In what ways is the area under aquaculture related to salinity prevalence and history of storm surges,
# accounting for distance and elevation from the sea. 

# Model 1: Basic pooled OLS regression
model_1 <- lm(Aqua_ha ~ postSurge + Saline_perc + 
               NEAR_DIST + DEM_avg, data = data)
summary(model_1)

# Model 2: Adding the interaction term
model_2 <- lm(Aqua_ha ~ postSurge + Saline_perc + 
               salinity_storm_interaction + NEAR_DIST + 
               DEM_avg, data = data)
summary(model_2)


# Model 3 - State level FE + Two way clustering (also clustering for villages across years)
model_3 <- feols(Aqua_ha ~ postSurge + Saline_perc + 
                   salinity_storm_interaction + 
                   DEM_avg | State + Sea_Dist, 
                 data = data, 
                 vcov = ~UniqueID + District)  # Two-way clustering
summary(model_3)


models <- list()
models[['Basic OLS']] <- model_1
models[['OLS + Interaction']] <- model_2
models[['Fixed Effects (State+Dist)']] <- model_3
msummary(models, stars = c('*' = .1, '**' = .05, '***' = .01),gof_omit=c("BIC|AIC|RMSE|R2 Within Adj."),coef_omit=c("(Intercept)"), filename = 'table.rtf')


# For persistent Salinity - Using the 5-year average salinity
# Model 4: Basic pooled OLS regression
model_4 <- lm(Aqua_ha ~ postSurge + avg_salinity_5yr + 
                NEAR_DIST + DEM_avg, data = data)
summary(model_4)

# Model 5: Adding the interaction term
model_5 <- lm(Aqua_ha ~ postSurge + avg_salinity_5yr + 
                avg_salinity_5yr*postSurge + NEAR_DIST + 
                DEM_avg, data = data)
summary(model_5)


# Model 6: State level FE + Two way clustering (also clustering for villages across years)
model_6 <- feols(Aqua_ha ~ postSurge + avg_salinity_5yr + 
                   avg_salinity_5yr*postSurge + DEM_avg | State + Sea_Dist, data = data, vcov=~UniqueID + District)
summary(model_6)

models <- list()
models[['Basic OLS']] <- model_4
models[['OLS + Interaction']] <- model_5
models[['Fixed Effects (State+Dist)']] <- model_6
msummary(models, stars = c('*' = .1, '**' = .05, '***' = .01),gof_omit=c("BIC|AIC|RMSE|R2 Within Adj."),coef_omit=c("(Intercept)"), filename = 'table.rtf')





# With Salinity as binary 
# Model 7: Basic pooled OLS regression
model_7 <- lm(Aqua_ha ~ postSurge + Saline + 
                NEAR_DIST + DEM_avg, data = data)
summary(model_7)

# Model 8: Adding the interaction term
model_8 <- lm(Aqua_ha ~ postSurge + Saline + 
                Saline * postSurge + NEAR_DIST + 
                DEM_avg, data = data)
summary(model_8)


# Model 9 - State level FE + Two way clustering (also clustering for villages across years)
model_9 <- feols(Aqua_ha ~ postSurge + Saline + 
                   Saline * postSurge + 
                   DEM_avg | State + Sea_Dist, 
                 data = data, 
                 vcov = ~UniqueID + District)  # Two-way clustering
summary(model_9)

models <- list()
models[['Basic OLS']] <- model_7
models[['OLS + Interaction']] <- model_8
models[['Fixed Effects (State+Dist)']] <- model_9

msummary(models, stars = c('*' = .1, '**' = .05, '***' = .01),gof_omit=c("BIC|AIC|RMSE|R2 Within Adj."),coef_omit=c("(Intercept)"), filename = 'table.rtf')


# All FE (State and Distance) + Two-way clustered SEs
models <- list()
models[['Salinity Area']] <- model_3
models[['Persistent Salinity']] <- model_6
models[['Above median Saline']] <- model_9

msummary(models, stars = c('*' = .1, '**' = .05, '***' = .01),gof_omit=c("BIC|AIC|RMSE|R2 Within Adj."),coef_omit=c("(Intercept)"), filename = 'table.rtf')



# Two-way fixed effects 
# Note: Two-way fixed effect controls for village-level fixed effects (absorbing time-invariant differences between villages) - which is why NEAR_DIST and DEM_avg get dropped.
# Since we still don't have pre-storm years for any storm affected villages, this postStorm is also getting dropped. 
# It controls for year fixed effects (accounting for changes common to all villages in a given year).
# The model focuses on within-village variation over time rather than cross-sectional differences.
# ! Do this AFTER adding Landsat5 years 

# Model 5: Fixed effects model (two-way fixed effect - village and year fixed effects)
model_5 <- plm(Aqua_ha ~ postSurge + Saline_perc + 
                 postSurge * Saline_perc, 
               data = pdata, 
               model = "within", 
               effect = "twoways")
summary(model_5)

# Testing for serial correlation
pbgtest(model_5)

# Very high test statistic and very low p-value suggests that we reject the null hypothesis of no serial correlation 
# Serial correlation violates a key assumption of standard regression models
# standard errors are likely biased (usually underestimated)
# significance tests may be invalid (i.e. "significant" effects that aren't really significant)

# Since serial correlation is present, we could: 
# 1. Use robust standard errors with vcovHC()
# 2. Consider including lagged dependent variables in the model 
# 3. Consider using first differences or other dynamic panel data methods
# 4. Try Newey-West standard errors or other corrections

# Using robust standard errors
robust_model_5 <- coeftest(model_5, vcov = vcovHC(model_5, type = "HC1"))
print(robust_model_5)





# ----------------------
# VISUALIZATIONS
# ----------------------

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
  Saline_perc = seq(min(data$Saline_perc, na.rm = TRUE),
                            max(data$Saline_perc, na.rm = TRUE),
                            length.out = 100),
  postSurge = c(0, 1),
  NEAR_DIST = mean(data$NEAR_DIST, na.rm = TRUE),
  DEM_avg= mean(data$DEM_avg, na.rm = TRUE)
)

grid_data$salinity_storm_interaction <- grid_data$Saline_perc * grid_data$postSurge
grid_data$predicted <- predict(model_2, newdata = grid_data)

# Plot the interaction effect
ggplot(grid_data, aes(x = Saline_perc, y = predicted, color = factor(postSurge))) +
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
