
###########################################
# CoSal for the Geomedian Composite Data #
###########################################

library(caret)
library(caretEnsemble)
library(pROC)
library(data.table)
library(dplyr)
library(scales)
library(glmnet)
library(xgboost)
library(tidyverse)
library(ggplot2)
library(gtsummary) # For summary tables
library(modelsummary)# For summary tables
library(mgcv) # GAM model fit 
library(randomForest) # to apply machine learning frameworks
library(datawizard) # for normalize()
library(nnet) # For ANN
library(neuralnet) # For more control on the architecture of ANN
library(glmnet) # For lasso regression 
library(MASS) # for stepwise regression
library(scales) # for scaling data from 0 to 1
library(kernlab)
library(naivebayes)

# Set working directory 
setwd("C:/Users/garim/Documents/GitHub_LocalRepository/aquaculture_patterns")
getwd()


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Step 1. Read the csv with soil field data and band information (extracted from ArcGIS)
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
soil_data <- read.csv("data/soil_data_raw.csv")

soil_data$Blue <- soil_data$b1_Band
soil_data$Green <- soil_data$b2_Band
soil_data$Red <- soil_data$b3_Band
soil_data$NIR <- soil_data$b4_Band
soil_data$SWIR1 <- soil_data$b5_Band
soil_data$SWIR2 <- soil_data$b6_Band
soil_data$NDVI <- soil_data$b7_Band
soil_data$NDWI <- soil_data$b8_Band
soil_data$NDSI <- soil_data$b9_Band

# Remove observations with missing Band values from Soil Data (9 points are on the masked areas with the geomedian classified image vs 7 in percentile composite)
soil_data <- soil_data[!is.na(soil_data$Blue), ]

# Keep soil moisture data for reflectance assessment
head(soil_data)
summary(soil_data$Soil_Moisture_Lab) # have field measures to fill the 42 missing lab measures, but leaving out for now. 

# Make a clean subset of data
soil_data <- soil_data[, c("Name", "EC", "pH", "Soil_Moisture_Lab", "Blue", "Green", "Red", "NIR", "SWIR1", "SWIR2")]
head(soil_data)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Step 2. Create Indices 
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# Convert EC into binary (high low salinity - threshold at 1900) 
soil_data$EC_all <- soil_data$EC # create a copy of continuous EC values - EC will take on different values for analysis 
soil_data$EC_bin <- ifelse(soil_data$EC >= 1900, 1, 0)
table(soil_data$EC_bin)

# Create Normalised differences for all bands apart from the regularly known indices
# Blue and red 
soil_data$NBR <- (soil_data$Blue-soil_data$Red)/(soil_data$Blue+soil_data$Red)

# Blue and green 
soil_data$NBG <- (soil_data$Blue-soil_data$Green)/(soil_data$Blue+soil_data$Green)

# Blue and NIR
soil_data$NBNIR <- (soil_data$Blue-soil_data$NIR)/(soil_data$Blue+soil_data$NIR)

# Blue and SWIR1
soil_data$NBSWIR1 <- (soil_data$Blue-soil_data$SWIR1)/(soil_data$Blue+soil_data$SWIR1)

# Blue and SWIR2 
soil_data$NBSWIR2 <- (soil_data$Blue-soil_data$SWIR2)/(soil_data$Blue+soil_data$SWIR2)

# Red and Green (also NDVI) 
soil_data$NDVI <- (soil_data$Red-soil_data$Green)/(soil_data$Red+soil_data$Green)

# Red and NIR (NDSI2 or Normalised Difference Salinity Index 2 as per as per Khan et al 2001 in Nguyen et al 2020)
soil_data$NDSI2 <- (soil_data$Red-soil_data$NIR)/(soil_data$Red+soil_data$NIR)

# Red and SWIR1
soil_data$NRSWIR1 <- (soil_data$Red-soil_data$SWIR1)/(soil_data$Red+soil_data$SWIR1)

# Red and SWIR2
soil_data$NRSWIR2 <- (soil_data$Red-soil_data$SWIR2)/(soil_data$Red+soil_data$SWIR2)

# Green and NIR (also NDWI) 
soil_data$NDWI <- (soil_data$Green-soil_data$NIR)/(soil_data$Green+soil_data$NIR)

# Green and SWIR1
soil_data$NGSWIR1 <- (soil_data$Green-soil_data$SWIR1)/(soil_data$Green+soil_data$SWIR1)

# Green and SWIR2
soil_data$NGSWIR2 <- (soil_data$Green-soil_data$SWIR2)/(soil_data$Green+soil_data$SWIR2)

# NIR and SWIR1
soil_data$NNIRSWIR1 <- (soil_data$NIR-soil_data$SWIR1)/(soil_data$NIR+soil_data$SWIR1)

# NIR and SWIR2 
soil_data$NNIRSWIR2 <- (soil_data$NIR-soil_data$SWIR2)/(soil_data$NIR+soil_data$SWIR2)

# SWIR1 and SWIR2 (also NDSI as per the Index Database: https://www.indexdatabase.de/db/is.php?sensor_id=168 )
soil_data$NDSI1 <- (soil_data$SWIR1-soil_data$SWIR2)/(soil_data$SWIR1+soil_data$SWIR2)



# Calculate other widely used Salinity indices (listed in Nguyen et al. 2020)
# Salinity Index 1 = sqrt(green^2+red^2)
soil_data$SI1 <- sqrt((soil_data$Green)^2 + (soil_data$Red)^2) 

# Salinity Index 2 = sqrt(green x red)
soil_data$SI2 <- sqrt(soil_data$Green * soil_data$Red)

# Salinity Index 3 = sqrt(blue x red) 
soil_data$SI3 <- sqrt(soil_data$Blue * soil_data$Red)

# salinity index 4 = red x NIR / green 
soil_data$SI4 <- (soil_data$Red * soil_data$NIR / soil_data$Green)  

# salinity index 5 = blue/red 
soil_data$SI5 <- (soil_data$Blue / soil_data$Red)   

# Soil Adjusted Vegetation Index (SAVI) = ((1.5)x NIR) - (red/0.5) + NIR + Red 
soil_data$SAVI <- (1.5 * soil_data$NIR) - (0.5 * soil_data$Red) + soil_data$NIR + soil_data$Red

# Vegetation Soil Salinity Index (VSSI) = (2 x green) - 5 x (red + NIR) 
soil_data$VSSI <- (2 * soil_data$Green) - 5 * (soil_data$Red + soil_data$NIR)





# Keep only necessary columns
soil_data <- soil_data[, c("Name", "EC", "EC_bin", "pH", "Soil_Moisture_Lab", "Blue", "Green", "Red", "NIR", "SWIR1", "SWIR2",
                           "NDVI", "NDWI", "NDSI1", "NDSI2", "SI1", "SI2", "SI3", "SI4", "SI5", "SAVI", "VSSI",
                           "NBR", "NBG", "NBNIR", "NBSWIR1", "NBSWIR2", "NRSWIR1", "NRSWIR2", "NGSWIR1", "NGSWIR2",
                           "NNIRSWIR1", "NNIRSWIR2")]

# Keep only numeric fields
soil_data_numeric <- soil_data[, sapply(soil_data, is.numeric)]

head(soil_data)

write.csv(soil_data, "data/soil_data_allindices.csv")


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Step 3. CoSal Model
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# ================ 1. Load functions ===============
# Load accuracy metrics function
source("codes/0_AccuracyFunction.R")


# ================ 2. Reduce data dimensionality  ===============

# Lasso variable selection 
set.seed(123)
drop.cols <- c('EC')
x <- soil_data_numeric %>% dplyr::select(-one_of(drop.cols))
# scale x
for(i in 1:ncol(x)){
  x[, i] <- rescale(x[, i])
}

x <- data.matrix(x)

y <- soil_data_numeric$EC
table(y)

cvmodel <- cv.glmnet(x, y, alpha=1, family='binomial') # did not converge
png(filename='outputs/lasso_lambda.png', width=800, height=400, res=96)
plot(cvmodel)
dev.off()

best_lambda <- cvmodel$lambda.min
best_lambda

# we can also tweak lambda to see
bestlasso <- glmnet(x, y, alpha=1, lambda=best_lambda, family='binomial')
coef(bestlasso)
# ! Update selected predictors
# With percentile composite: Green, SWIR2, NDWI, NDSI1, VSSI, NBNIR, NBSWIR2, NRSWIR1, NGSWIR1, NNIRSWIR1
# With geomedian composite: Red, SWIR2, SI4, VSSI, NBNIR, NBSWIR2, NRSWIR2, NGSWIR1

soil_data_numeric <- dplyr::select(soil_data_numeric,
                                   c('EC', 'Red', 'SWIR2', 'SI4' ,'VSSI', 'NBNIR', 'NBSWIR2', 'NRSWIR2', 'NGSWIR1'))


# ================ 3. Stratified Train-Test Split ===============
# Convert EC to factor for model classification
soil_data_numeric$EC <- as.factor(soil_data_numeric$EC)


set.seed(123)
train_indices <- createDataPartition(soil_data_numeric$EC, p = 0.8, list = FALSE)
train_data <- soil_data_numeric[train_indices, ]
test_data <- soil_data_numeric[-train_indices, ]

# Convert EC (0,1) to factor levels X0 and X1
levels(train_data$EC) <- make.names(levels(train_data$EC))
levels(test_data$EC) <- make.names(levels(test_data$EC))
train_data$EC <- factor(train_data$EC, levels = c("X0", "X1"))
test_data$EC <- factor(test_data$EC, levels = c("X0", "X1"))

# Debug: Check class balance
print("Train class distribution:")
print(table(train_data$EC))
print("Test class distribution:")
print(table(test_data$EC))

# ================ 4. Remove Near-Zero Variance Predictors ===============
nzv <- nearZeroVar(train_data, saveMetrics = TRUE)
if (any(nzv$nzv)) {
  print("Near-zero variance predictors found and removed:")
  print(nzv[nzv$nzv, ])
  train_data <- train_data[, !nzv$nzv]
  test_data <- test_data[, !nzv$nzv]
}

# Define cross-validation control
train_control <- trainControl(method = "cv", number = 5, savePredictions = "final", classProbs = TRUE)
# Using cv with 5 fold here. Could also test method = "LOOCV" (Leave-One-Out Cross-Validation) 


# ================ 5. Train Models Individually ===============
li_models <- c("rf", "rpart", "nnet", "svmRadial", "gbm", "naive_bayes", "xgbTree", "knn", "glmnet")

li_multi <- caretList(
  EC ~ .,
  data = train_data,
  trControl = train_control,
  methodList = li_models
)

# Warning comes from xgboost

# Check models - If null remove 
print(li_multi)



# ================ 6. Model Performance & Selection ===============
model_performances <- resamples(li_multi)
results <- summary(model_performances)
print(results)

# Extract accuracies
accuracies <- unlist(lapply(li_multi, function(model) {
  max(model$results$Accuracy)
}))

# Assess detailed statistics
performance_stats <- summary(model_performances)$statistics

# Use statistical methods to set a threshold
# Method 1: Mean threshold
mean_threshold <- mean(accuracies)

# # Method 2: Median threshold
# median_threshold <- median(accuracies)
# 
# # Method 3: Mean minus one SD
# sd_threshold <- mean(accuracies) - sd(accuracies)

selected_models <- names(accuracies[accuracies > mean_threshold])

print(paste("Mean accuracy:", round(mean_threshold, 4)))
print(paste("Models above threshold:", paste(selected_models, collapse=", ")))


# Visual inspection to find  a natural breakpoint?
# Plot model performances
model_comparison <- dotplot(model_performances)
plot(model_comparison)



# ================ 7. Create the Stacked Ensemble ===============
MultiStratEnsemble <- caretStack(
  li_multi,
  method = "rf",  
  metric = "ROC",
  trControl = trainControl(
    method = "cv", 
    number = 5, 
    classProbs = TRUE,
    summaryFunction = twoClassSummary
  )
)

# Check the ensemble object
print(summary(MultiStratEnsemble))



# ================ 8. Make Predictions & Fix Type Issues ===============
# Predict probabilities
train_preds <- predict(MultiStratEnsemble, newdata=train_data)
test_preds <- predict(MultiStratEnsemble, newdata=test_data)

# Ensure output is a data frame with X0 and X1 probabilities
print(str(train_preds))  
print(str(test_preds))  

# Extract probabilities of class X1
train_preds_numeric <- train_preds$X1
test_preds_numeric <- test_preds$X1

# Check lengths to avoid errors
print(length(train_preds_numeric))
print(length(test_preds_numeric))
print(length(train_data$EC))
print(length(test_data$EC))


# ================ 9. Determine Best Threshold (Youden's J) ===============
test_ec_numeric <- as.numeric(test_data$EC == "X1")
roc_obj <- roc(test_ec_numeric, test_preds_numeric)
best_threshold_df <- coords(roc_obj, "best", best.method = "youden")
best_threshold <- best_threshold_df$threshold

print(paste("Best threshold:", round(best_threshold, 4)))
# Best threshold: 0.475

# test sensitivity to thresholds
li_thresholds <- seq(0.1, 0.9, 0.01)
df <- data.frame(matrix(ncol=4, nrow=0))
colnames(df) <- c("threshold", "train_acc","test_acc", "diff")

for(threshold in li_thresholds){
  # Predict classes 
  train_predicted_class <- ifelse(train_preds[, "X1"] > threshold, "X1", "X0")
  test_predicted_class <- ifelse(test_preds[, "X1"] > threshold, "X1", "X0")
  
  # Calculate metrics
  train_ensemble_metrics <- calculate_classification_metrics(train_data$EC, train_predicted_class)
  train_acc <- train_ensemble_metrics$Accuracy
  test_ensemble_metrics <- calculate_classification_metrics(test_data$EC, test_predicted_class)
  test_acc <- test_ensemble_metrics$Accuracy
  diff <- train_acc-test_acc
  
  df[nrow(df)+1, ] <- c(threshold, train_acc, test_acc, abs(diff))
  
  #  if(abs(diff)<0.1){
  #    print(threshold)
  #  }
}

write.csv(df, 'supplementary/2_Threshold_tests.csv')

# min(df$diff)

mean(df$diff)

# ================ 10. Prediction ===============

# Convert probabilities to class predictions
train_predicted_class <- ifelse(train_preds_numeric > best_threshold, "X1", "X0")
test_predicted_class <- ifelse(test_preds_numeric > best_threshold, "X1", "X0")

# Ensure predictions match EC factor levels
train_predicted_class <- factor(train_predicted_class, levels = c("X0", "X1"))
test_predicted_class <- factor(test_predicted_class, levels = c("X0", "X1"))

# Final length check before evaluation
print(length(train_predicted_class))
print(length(test_predicted_class))



# ================ 11. Evaluate Model Performance ===============
train_ensemble_metrics <- calculate_classification_metrics(train_data$EC, train_predicted_class)
test_ensemble_metrics <- calculate_classification_metrics(test_data$EC, test_predicted_class)

print("Train Metrics:")
print(train_ensemble_metrics)

print("Test Metrics:")
print(test_ensemble_metrics)



# ================ 12. Save Model ===============
saveRDS(MultiStratEnsemble, "ensemble_geomedian.rds")
test <- readRDS("ensemble_geomedian.rds")


