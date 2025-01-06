# Required packages
library(caret)
library(glmnet)
library(xgboost)
library(randomForest)
library(nnet)
library(pROC)

#' First Layer Model Training Function
#' @param feature_data List of feature matrices 
#' @param labels Vector of binary labels
#' @param k Number of top models to select per feature
#' @return List of trained base models
train_first_layer <- function(feature_data, labels, k = 5) {
  base_models <- list()
  
  for(feature_name in names(feature_data)) {
    feature_matrix <- feature_data[[feature_name]]
    
    # Parameter grids for each algorithm
    param_grids <- list(
      glm = expand.grid(
        alpha = seq(0, 1, 0.2)
      ),
      xgboost = expand.grid(
        nrounds = c(100, 200),
        max_depth = c(3, 6),
        eta = c(0.01, 0.1),
        gamma = c(0, 0.1),
        subsample = c(0.7, 0.8, 0.9)
      ),
      rf = expand.grid(
        mtry = c(2, 4, 6),
        ntree = c(500, 1000)
      ),
      nn = expand.grid(
        size = c(3, 4, 5),
        decay = c(0.001, 0.01, 0.1)
      )
    )
    
    # Store results for each model
    results <- list()
    
    # GLM
    for(i in 1:nrow(param_grids$glm)) {
      set.seed(123)
      # Let cv.glmnet choose lambda sequence automatically
      cv_fit <- cv.glmnet(as.matrix(feature_matrix), labels,
                         alpha = param_grids$glm$alpha[i],
                         family = "binomial")
      # Use lambda.min for predictions
      pred <- predict(cv_fit, as.matrix(feature_matrix), 
                     s = "lambda.min", type = "response")
      auc <- auc(roc(labels, as.vector(pred)))
      results[[paste0("glm_", i)]] <- list(model = cv_fit, auc = auc, 
                                          params = param_grids$glm[i,])
    }
    
    # XGBoost
    for(i in 1:nrow(param_grids$xgboost)) {
      set.seed(123)
      dtrain <- xgb.DMatrix(data = as.matrix(feature_matrix), label = labels)
      model <- xgboost(data = dtrain,
                      nrounds = param_grids$xgboost$nrounds[i],
                      max_depth = param_grids$xgboost$max_depth[i],
                      eta = param_grids$xgboost$eta[i],
                      gamma = param_grids$xgboost$gamma[i],
                      subsample = param_grids$xgboost$subsample[i],
                      objective = "binary:logistic",
                      verbose = 0)
      pred <- predict(model, as.matrix(feature_matrix))
      auc <- auc(roc(labels, pred))
      results[[paste0("xgb_", i)]] <- list(model = model, auc = auc, 
                                          params = param_grids$xgboost[i,])
    }
    
    # Random Forest
    for(i in 1:nrow(param_grids$rf)) {
      set.seed(123)
      model <- randomForest(x = feature_matrix, y = as.factor(labels),
                          mtry = param_grids$rf$mtry[i],
                          ntree = param_grids$rf$ntree[i])
      pred <- predict(model, feature_matrix, type = "prob")[,2]
      auc <- auc(roc(labels, pred))
      results[[paste0("rf_", i)]] <- list(model = model, auc = auc, 
                                         params = param_grids$rf[i,])
    }
    
    # Neural Network
    for(i in 1:nrow(param_grids$nn)) {
      set.seed(123)
      model <- nnet(feature_matrix, labels,
                   size = param_grids$nn$size[i],
                   decay = param_grids$nn$decay[i],
                   maxit = 200,
                   trace = FALSE)
      pred <- predict(model, feature_matrix)
      auc <- auc(roc(labels, pred))
      results[[paste0("nn_", i)]] <- list(model = model, auc = auc, 
                                         params = param_grids$nn[i,])
    }
    
    # Select top k models for this feature type
    aucs <- sapply(results, function(x) x$auc)
    top_models <- results[order(aucs, decreasing = TRUE)[1:k]]
    base_models[[feature_name]] <- top_models
  }
  
  return(base_models)
}

#' Generate predictions from first layer models
#' @param base_models List of trained base models from first layer
#' @param feature_data List of feature matrices
#' @return Matrix of predictions from all base models
get_base_predictions <- function(base_models, feature_data) {
  n_samples <- nrow(feature_data[[1]])
  total_models <- sum(sapply(base_models, length))
  base_predictions <- matrix(0, nrow = n_samples, ncol = total_models)
  colnames(base_predictions) <- character(total_models)
  
  col_idx <- 1
  for(feature_name in names(base_models)) {
    feature_matrix <- feature_data[[feature_name]]
    feature_models <- base_models[[feature_name]]
    
    for(model_name in names(feature_models)) {
      model_info <- feature_models[[model_name]]
      model <- model_info$model

      if(inherits(model, "cv.glmnet")) {
        pred <- predict(model, newx = as.matrix(feature_matrix), 
                       s = "lambda.min", type = "response")
      } else if(inherits(model, "xgb.Booster")) {
        pred <- predict(model, as.matrix(feature_matrix))
      } else if(inherits(model, "randomForest")) {
        pred <- predict(model, feature_matrix, type = "prob")[,2]
      } else if(inherits(model, "nnet")) {
        pred <- predict(model, feature_matrix)
      }
      
      base_predictions[, col_idx] <- as.vector(pred)
      colnames(base_predictions)[col_idx] <- paste0(feature_name, "_", model_name)
      col_idx <- col_idx + 1
    }
  }
  
  return(base_predictions)
}

#' Second Layer Model Training Function
#' @param base_predictions Matrix of predictions from base models
#' @param labels Vector of binary labels
#' @param n_models Number of final models to ensemble
#' @return List of trained stacking models
train_second_layer <- function(base_predictions, labels, n_models = 10) {
  # Parameter grids for stacking algorithms
  param_grids <- list(
    glm = expand.grid(
      alpha = seq(0, 1, 0.2)
    ),
    xgboost = expand.grid(
      nrounds = c(100, 200),
      max_depth = c(3, 6),
      eta = c(0.01, 0.1)
    ),
    rf = expand.grid(
      mtry = c(2, 4, 6),
      ntree = c(500, 1000)
    )
  )
  
  # Store results for each stacking model
  stack_results <- list()
  
  # Train stacking models
  for(algorithm in c("glm", "xgboost", "rf")) {
    params <- param_grids[[algorithm]]
    
    for(i in 1:nrow(params)) {
      set.seed(123)
      
      if(algorithm == "glm") {
        model <- cv.glmnet(as.matrix(base_predictions), labels,
                          alpha = params$alpha[i],
                          family = "binomial")
        pred <- predict(model, as.matrix(base_predictions), 
                       s = "lambda.min", type = "response")
      } else if(algorithm == "xgboost") {
        dtrain <- xgb.DMatrix(data = as.matrix(base_predictions), label = labels)
        model <- xgboost(data = dtrain,
                        nrounds = params$nrounds[i],
                        max_depth = params$max_depth[i],
                        eta = params$eta[i],
                        objective = "binary:logistic",
                        verbose = 0)
        pred <- predict(model, as.matrix(base_predictions))
      } else if (algorithm == "nnet") {
			model <- nnet(feature_matrix, labels,
					   size = params$size[i],
					   decay = params$decay[i],
					   maxit = 200,
					   trace = FALSE)
		  pred <- predict(model, feature_matrix)
	  } else {
        model <- randomForest(x = base_predictions, y = as.factor(labels),
                            mtry = params$mtry[i],
                            ntree = params$ntree[i])
        pred <- predict(model, base_predictions, type = "prob")[,2]
      }
      
      # Calculate AUC
      auc <- auc(roc(labels, as.vector(pred)))
      stack_results[[paste0(algorithm, "_", i)]] <- list(
        model = model,
        auc = auc,
        params = params[i,]
      )
    }
  }
  
  # Select top n_models
  aucs <- sapply(stack_results, function(x) x$auc)
  top_models <- stack_results[order(aucs, decreasing = TRUE)[1:n_models]]
  
  return(top_models)
}
