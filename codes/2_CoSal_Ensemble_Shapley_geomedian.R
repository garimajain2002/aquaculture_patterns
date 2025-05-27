# ______________________________________________________________________________
# ~~ Stack Ensemble multiple models with Binary EC and stratified sampling ~~ #
# ______________________________________________________________________________

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
library(caret)  # For bagging
library(MASS) # for stepwise regression
library(dplyr)
library(scales) # for scaling data from 0 to 1
library(caretEnsemble) # For ensemble
library(kernlab)
library(naivebayes)
library(pROC)

# ================ 1. Read data ===============

setwd("C:/Users/garim/Documents/GitHub_LocalRepository/aquaculture_patterns")
getwd()

source("codes/0_AccuracyFunction.R")

MultiStratEnsemble <- readRDS("ensemble_geomedian.rds") 

soil_data <- read.csv("data/soil_data_allindices.csv") # make sure soil_data has EC values

head(soil_data)

soil_data_numeric$EC <- as.factor(soil_data_numeric$EC_bin)

table(soil_data_numeric$EC)


# ================ 2. Prepare data ===============

# Select the type of EC to be used here
head(soil_data)

# Clean dataset
soil_data <- soil_data[, c("Name", "EC", "Blue", "Green", "Red", "NIR", "SWIR1", "SWIR2",
                           "NDVI", "NDWI", "NDSI1", "NDSI2", "SI1", "SI2", "SI3", "SI4", "SI5", "SAVI", "VSSI",
                           "NBR", "NBG", "NBNIR", "NBSWIR1", "NBSWIR2", "NRSWIR1", "NRSWIR2", "NGSWIR1", "NGSWIR2",
                           "NNIRSWIR1", "NNIRSWIR2")]


# Prepare soil data with only numeric fields
soil_data_numeric <- soil_data[, sapply(soil_data, is.numeric)]
soil_data_numeric$EC <- as.factor(soil_data_numeric$EC)

table(soil_data_numeric$EC)


# ================ 3. Stratified Sampling ===============

#### With Stratified sampling for training and testing data along with 5-fold###

  # Stratified sampling for train-test split
  set.seed(123)
  train_indices <- createDataPartition(soil_data_numeric$EC, p = 0.8, list = FALSE)
  train_data <- soil_data_numeric[train_indices, ]
  test_data <- soil_data_numeric[-train_indices, ]

  
#   # Convert EC into Factor (0 and 1 to X0 and X1)
#   levels(train_data$EC) <- make.names(levels(train_data$EC))
#   levels(test_data$EC) <- make.names(levels(test_data$EC))
#   train_data$EC <- factor(train_data$EC, levels = c("X0", "X1"))
#   test_data$EC <- factor(test_data$EC, levels = c("X0", "X1"))
#   
#   # Debug: Check class balance
#   print(paste("Train class distribution"))
#   print(table(train_data$EC))
#   print(paste("Test class distribution"))
#   print(table(test_data$EC))
#   
#   # Check and remove near-zero variance predictors
#   nzv <- nearZeroVar(train_data, saveMetrics = TRUE)
#   if (any(nzv$nzv)) {
#     print("Near-zero variance predictors found:")
#     print(nzv[nzv$nzv, ])
#     train_data <- train_data[, !nzv$nzv]
#     test_data <- test_data[, !nzv$nzv]
#   }
#   
#   # Define training control
#   train_control <- trainControl(method = "cv", number = 5, savePredictions = "final", classProbs = TRUE)
#   # Using cv with 5 fold here. Could also test method = "LOOCV" (Leave-One-Out Cross-Validation) 
#   
#   
#   # ================ 4. Train models ===============
#   
#   # Train multiple models using caretList
#   li_multi <- caretList(
#     EC ~ Blue_R + Red_R + Green_R + NIR_R + SWIR1_R + SWIR2_R + 
#       NDVI + NDWI + NDSI1 + NDSI2 + SI1 + SI2 + SI3 + SI4 + SI5 + SAVI + VSSI + 
#       NBR + NBG + NBNIR + NBSWIR1 + NBSWIR2 + NRSWIR1 + NRSWIR2 + NGSWIR1 + NGSWIR2 + NNIRSWIR1 + NNIRSWIR2,
#     data = train_data,
#     trControl = train_control,
#     methodList = c("rf", "rpart", "nnet", "svmRadial", "gbm", "naive_bayes", "xgbTree", "knn", "glmnet")  
#   )
#   
#   # Check the list of models that they are functioning and not null. If null remove. 
#   li_multi
#   
#   # Model/Feature selection
#   # Get performance metrics for all models
#   model_performances <- resamples(li_multi)
#   results <- summary(model_performances)
#   
#   print(results)
#   
#   accuracies <- unlist(lapply(li_multi, function(model) {
#     max(model$results$Accuracy)
#   }))
#   
#   # Assess detailed statistics
#   performance_stats <- summary(model_performances)$statistics
#   
#   # Use statistical methods to set a threshold
#   # Method 1: Mean threshold
#   mean_threshold <- mean(accuracies)
#   
#   # Method 2: Median threshold
#   median_threshold <- median(accuracies)
#   
#   # Method 3: Mean minus one SD
#   sd_threshold <- mean(accuracies) - sd(accuracies)
#   
#   selected_models <- names(accuracies[accuracies > mean_threshold])
#   
#   print(paste("Mean accuracy:", round(mean_threshold, 4)))
#   print(paste("Models above threshold:", paste(selected_models, collapse=", ")))
#   
#   
#   # Visual inspection to find  a natural breakpoint?
#   # Plot model performances
#   model_comparison <- dotplot(model_performances)  
#   plot(model_comparison )
#   
#   
#   # ================ 5. Create the stack ensemble ===============
#   # Define stack ensemble
#   MultiStratEnsemble <- caretStack(
#     li_multi,
#     method = "glm",  # Stacking method
#     metric = "Accuracy",
#     trControl = trainControl(method = "cv", number = 5, classProbs = TRUE)
#   )
#   
#   # # Check the ensemble object
#   # MultiStratEnsemble
#   # #str(MultiStratEnsemble)
#   
#   
#   # ================ 6. Evaluate the Stack Ensemble ===============
#   #   # Summary of the ensemble model
#   # summary(MultiStratEnsemble)
#   
#   # Predict probabilities
#   # Check individual model predictions ("rf", "rpart", "nnet", "svmRadial", "gbm", "naive_bayes", "xgbTree", "knn", "glmnet")
#   rf_prob <- predict(li_multi$rf, newdata = test_data, type = "prob")
#   class(li_multi$rf)
#   rpart_prob <- predict(li_multi$rpart, newdata = test_data, type = "prob")
#   nnet_prob <- predict(li_multi$nnet, newdata = test_data, type = "prob")
#   svmRadial_prob <- predict(li_multi$svmRadial, newdata = test_data, type = "prob")
#   gbm_prob <- predict(li_multi$gbm, newdata = test_data, type = "prob")
#   naive_bayes_prob <- predict(li_multi$naive_bayes, newdata = test_data, type = "prob")
#   xgb_prob <- predict(li_multi$xgbTree, newdata = test_data, type = "prob")
#   knn_prob <- predict(li_multi$knn, newdata = test_data, type = "prob")
#   glmnet_prob <- predict(li_multi$glmnet, newdata = test_data, type = "prob")
#  
#   # Combine probabilities (average)
#   ensemble_probabilities <- (rf_prob + rpart_prob + nnet_prob + svmRadial_prob + gbm_prob + naive_bayes_prob + xgb_prob + knn_prob + glmnet_prob) / 9
#   
#   
#   # determine the best threshold
#   # Create ROC object
#   roc_obj <- roc(test_data$EC, ensemble_probabilities[, "X1"])
#   
#   # Find best threshold using Youden's J statistic
#   best_threshold <- coords(roc_obj, "best", best.method = "youden")
#   
#   # Alternative: find threshold that maximizes specificity + sensitivity
#   best_threshold_alt <- coords(roc_obj, "best", best.method = "closest.topleft")
#   
#   # Print the threshold
#   print(best_threshold)
#  
#   
#     # Predict classes 
#   predicted_class <- ifelse(ensemble_probabilities[, "X1"] > best_threshold$threshold, 1, 0)
#   
#   
#   # Calculate metrics
#   ensemble_metrics <- calculate_classification_metrics(test_data$EC, predicted_class)
#   
#   
#   # Print ensemble metrics
#   print(ensemble_metrics)
  
  
  # -----------------------------------------------------------
  # Calculate and visualise Shapley values for all parameters
  # ------------------------------------------------------------
  
  library(iml)
  library(ggplot2)
  library(tidyverse)
  
  # Define the custom prediction wrapper
  predict_wrapper <- function(model, newdata) {
    preds <- predict(model, newdata = newdata)
    # Convert to numeric matrix with class probabilities
    return(as.matrix(preds))
  }
  
  # Function to calculate Shapley values
  calculate_shapley <- function(model, train_data, test_data, n_samples = 50) {
    # Create predictor object
    predictor <- Predictor$new(
      model = model,
      data = train_data[, -which(names(train_data) == "EC")],  # Remove target variable
      y = train_data$EC,
      predict.fun = predict_wrapper
    )
    
    results <- list()
    feature_names <- setdiff(names(train_data), "EC")
    
    # Calculate Shapley values for each test instance
    for(i in 1:nrow(test_data)) {
      cat(sprintf("Computing Shapley values for instance %d of %d\n", i, nrow(test_data)))
      
      shapley <- Shapley$new(
        predictor = predictor,
        x.interest = test_data[i, feature_names],
        sample.size = n_samples
      )
      
      results[[i]] <- shapley$results
    }
    
    return(results)
  }
  
  # Function to aggregate and analyze Shapley values
  analyze_shapley_results <- function(shapley_results) {
    # Combine all results into a single dataframe
    all_results <- do.call(rbind, shapley_results)
    
    # Calculate average absolute Shapley values for each feature
    feature_importance <- all_results %>%
      group_by(feature) %>%
      summarise(
        mean_abs_effect = mean(abs(phi)),
        sd_effect = sd(phi)
      ) %>%
      arrange(desc(mean_abs_effect))
    
    return(feature_importance)
  }
  
  # Function to plot Shapley results
  plot_shapley_importance <- function(feature_importance) {
    ggplot(feature_importance, aes(x = reorder(feature, mean_abs_effect), y = mean_abs_effect)) +
      geom_bar(stat = "identity", fill = "steelblue") +
      geom_errorbar(aes(ymin = mean_abs_effect - sd_effect, 
                        ymax = mean_abs_effect + sd_effect), 
                    width = 0.2) +
      coord_flip() +
      theme_minimal() +
      labs(
        title = "Feature Importance Based on Shapley Values",
        x = "Features",
        y = "Mean Absolute Shapley Value"
      )
  }
  
  # Execute the analysis
  # 1. Calculate Shapley values
  shapley_results <- calculate_shapley(MultiStratEnsemble, train_data, test_data, n_samples = 50)
  
  # 2. Analyze results
  feature_importance <- analyze_shapley_results(shapley_results)
  print(feature_importance)
  
  # 3. Plot results
  importance_plot <- plot_shapley_importance(feature_importance)
  print(importance_plot)
  
  # 4. Save results if needed
  write.csv(feature_importance, "outputs/shapley_feature_importance_ensemble_lassoPredictors.csv", row.names = FALSE)
  ggsave("outputs/shapley_importance_plot.png", importance_plot, width = 10, height = 8)
  
  # If you want to examine individual predictions:
  example_instance <- 1  # Change this to look at different test instances
  print("Shapley values for a single prediction:")
  print(shapley_results[[example_instance]] %>% arrange(desc(abs(phi))))

    