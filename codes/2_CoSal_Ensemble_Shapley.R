# -----------------------------------------------------------
# Calculate and visualise Shapley values for all parameters
# ------------------------------------------------------------
# # ~~ Stack Ensemble multiple models with Binary EC and stratified sampling ~~ #
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
library(iml)


# ================ 1. Read data and models ===============
setwd("C:/Users/garim/Documents/GitHub_LocalRepository/aquaculture_patterns")
getwd()

source("codes/0_AccuracyFunction.R")

MultiStratEnsemble <- readRDS("ensemble_geomedian.rds") 

soil_data <- read.csv("data/soil_data_allindices.csv")

head(soil_data)

soil_data_numeric$EC <- as.factor(soil_data_numeric$EC)

table(soil_data_numeric$EC)


# ================ 2. Shapley ===============

#### With Stratified sampling for training and testing data along with 5-fold###

  # Stratified sampling for train-test split
  set.seed(123)
  train_indices <- createDataPartition(soil_data_numeric$EC, p = 0.8, list = FALSE)
  train_data <- soil_data_numeric[train_indices, ]
  test_data <- soil_data_numeric[-train_indices, ]


  
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

    