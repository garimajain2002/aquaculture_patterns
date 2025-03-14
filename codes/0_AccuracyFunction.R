# Script to store frequently used functions
getwd()

# Define functions
calculate_classification_metrics <- function(actual, predicted) {
  # Ensure both are factors with same levels
  actual <- factor(actual, levels = c("X0", "X1"))
  predicted <- factor(predicted, levels = c("X0", "X1"))
  
  # Create confusion matrix
  confusion <- table(actual, predicted)
  
  # Extract values safely
  TN <- confusion["X0", "X0"]
  FP <- confusion["X0", "X1"]
  FN <- confusion["X1", "X0"]
  TP <- confusion["X1", "X1"]
  
  # Replace NA with 0
  TN <- ifelse(is.na(TN), 0, TN)
  FP <- ifelse(is.na(FP), 0, FP)
  FN <- ifelse(is.na(FN), 0, FN)
  TP <- ifelse(is.na(TP), 0, TP)
  
  # Calculate metrics
  accuracy <- (TP + TN) / (TN + FP + FN + TP)
  precision <- TP / (TP + FP)
  recall <- TP / (TP + FN)
  f1_score <- 2 * (precision * recall) / (precision + recall)
  
  return(list(
    Accuracy = accuracy,
    Precision = precision,
    Recall = recall,
    F1_Score = f1_score
  ))
}