# Author: Camila Alvarez-Silva

library(dplyr)
library(xgboost)
library(caret)
library(pROC)
library(Metrics)
library(reshape2)
library(ggplot2)
library(groupdata2)
library(xpectr)

# =============================================================================
# XGBoost Helper Functions
# -----------------------------------------------------------------------------
# Two main functions:
#   1. xgboost_train()   – trains XGBoost models per omic layer using repeated
#                        stratified 5-fold CV; returns predictions, AUC,
#                        feature importance, and top-N feature sub-models.
# =============================================================================


# =============================================================================
# 1. xgboost_train
# -----------------------------------------------------------------------------
# Args:
#   omics      – character vector of omic layer names to iterate over
#   alpha      – L1 regularisation per omic (parallel to `omics`)
#   feature    – data.frame with target label column named "feature"
#   seed       – random seed for reproducibility
#   top        – max number of top features to evaluate (2:top)
#   df.matrix  – full feature matrix (samples × features)
#   nfolds:(Integer)	Number of folds for cross-validation.
#   fecal.features: (Character Vector) List of microbiome-related features.
#   plasma.features: (Character Vector) 	List of host-related features.
#   metadata: (DataFrame)	Metadata containing sample IDs and clinical details.

# Returns: named list with models, predictions, AUC, importance, folds, etc.
# =============================================================================

xgboost_train <- function(omics, 
                          feature, 
                          seed = seed, 
                          top = top, 
                          df.matrix = df.matrix, 
                          nfolds = nfolds,
                          ntimes=ntimes, 
                          Microbiome.features = Microbiome.features, 
                          Host.features = Host.features, 
                          metadata = metadata) {
  {
  
  # Check if essential inputs exist
  if (missing(df.matrix) || missing(metadata) || missing(feature) || missing(nfolds) || missing(Microbiome.features) ||
      missing(Host.features) || missing(panel)) {
    stop("Missing required input")
  }
  
  # Shared XGBoost hyperparameters – edit here to affect all model calls
  xgb_params <- list(
    nround               = 100,
    alpha                = 15,   
    lambda               = 0,
    eta                  = 0.3,
    gamma                = 0,
    max_depth            = 6,
    early_stopping_rounds = 100,
    eval_metric          = "auc",
    objective            = "binary:logistic",
    nthread              = 20
  )
  
  # ---------------------------------------------------------------------------
  # Standard clinical/non-omics variables treated as single-feature models
  # ---------------------------------------------------------------------------
  standards.variables <- c(
    "standards_swe", "standards_te",   "standards_elf",  "standards_apri",
    "standards_at",  "standards_forns","standards_fib4", "standards_ft",
    "standards_cap", "standards_gpr",  "standards_aar",  "standards_proc3",
    "standards_apoa1","standards_nfs", "standards_cp",   "standards_meld",
    "standards_p3np","standards_ha"
  )
  
  # ---------------------------------------------------------------------------
  # Initialise output lists (one entry per omic)
  # ---------------------------------------------------------------------------
  loss.list                  <- list()
  cvAUC.list                 <- list()
  cvAUC.list.top             <- list()
  cvAUC.value                <- list()
  cvAUC.se                   <- list()
  folds.list                 <- list()
  score_list.all             <- list()
  predictions.list.all       <- list()
  predictions.list.all.top   <- list()
  model_list.all             <- list()
  model_list.all.top         <- list()
  feature.importance.list.all     <- list()
  feature.importance.list.top.all <- list()
  feature.importance.list2        <- list()   
  confusionMatrix.list.all        <- list()
  confusionMatrix.lit.top.all     <- list()
  out.list                   <- list()
  out.list.top               <- list()
  train.labels.list.out      <- list()
  val.labels.list.out        <- list()
  
  # ===========================================================================
  # MAIN LOOP – iterate over each omic layer
  # ===========================================================================
  for (j in seq_along(omics)) {
    
    omi        <- omics[j]
    # FIX: honour the per-omic alpha argument instead of hardcoding 15
    alpha.omic <- 15
    
    cat("####################################\n")
    cat(omi, "\n")
    cat("####################################\n")
    
    # -------------------------------------------------------------------------
    # Per-omic initialisation
    # -------------------------------------------------------------------------
    labelsValidation.list       <- list()
    score_list_inner            <- list()   
    predictions.list            <- list()
    predictions.list.top        <- list()
    model_list                  <- list()
    model_list.top              <- list()
    feature.importance.list     <- list()
    feature.importance.list.top <- list()
    confusionMatrix.list        <- list()
    confusionMatrix.list.top    <- list()
    train.labels.list           <- list()
    val.labels.list             <- list()
    
    # =========================================================================
    # STEP 1 – Subset feature matrix for the current omic layer
    # =========================================================================
    
    if (!omi %in% c("all", "Microbiome", "Host")) {
      omics.features  <- colnames(df.matrix)[startsWith(colnames(df.matrix), omi)]
      df.matrix.omics <- df.matrix[, omics.features, drop = FALSE]
      
    } else if (omi == "all") {
      df.matrix.omics <- df.matrix[, c(
        which(colnames(df.matrix) %in% Microbiome.features),
        which(colnames(df.matrix) %in% Host.features)
      )]
      
    } else if (omi == "Microbiome") {
      df.matrix.omics <- df.matrix[, colnames(df.matrix) %in% Microbiome.features,
                                   drop = FALSE]
      
    } else if (omi == "Host") {
      df.matrix.omics <- df.matrix[, colnames(df.matrix) %in% Host.features,
                                   drop = FALSE]
    }
    
    # Remove samples that are entirely NA
    if (omi %in% standards.variables) {
      df.matrix.omics <- df.matrix.omics[!is.na(df.matrix.omics[, omi]), ,
                                         drop = FALSE]
    } else {
      all_na_rows     <- apply(df.matrix.omics, 1, function(x) all(is.na(x)))
      df.matrix.omics <- df.matrix.omics[!all_na_rows, , drop = FALSE]
    }
    
    # Attach target label; use Row.names as row identifier
    df.matrix.omics            <- merge(feature, df.matrix.omics, by = 0)
    row.names(df.matrix.omics) <- df.matrix.omics$Row.names
    df.matrix.omics$Row.names  <- NULL
    
    train_data   <- as.matrix(df.matrix.omics[, colnames(df.matrix.omics) != "feature"])
    train_labels <- df.matrix.omics$feature
    
    # =========================================================================
    # STEP 2 – Create repeated stratified 5-fold CV indices (5 × 5 = 25 folds)
    # =========================================================================
    cat("STEP 2: Creating 5×5 stratified CV folds for", omi, "\n")
    
    nfolds <- nfolds
    ntimes <- ntimes
    folds  <- nkfold(y = train_labels, n = ntimes, k = nfolds,
                     stratified = TRUE, seed = seed, named = TRUE)
    
    folds.list[[omi]] <- folds
    
    # =========================================================================
    # STEP 3 – Tune on fold 1 and save loss curve for inspection
    # =========================================================================
    cat("STEP 3: Tuning on fold 1 and saving loss curve for", omi, "\n")
    
    x1 <- folds$Rep1Fold1
    dtrain.f1 <- xgb.DMatrix(
      data  = .fix_colnames(train_data[-x1, , drop = FALSE], omi, standards.variables),
      label = train_labels[-x1]
    )
    dval.f1 <- xgb.DMatrix(
      data  = .fix_colnames(train_data[x1,  , drop = FALSE], omi, standards.variables),
      label = train_labels[x1]
    )
    
    model.f1 <- xgb.train(
      data        = dtrain.f1,
      watchlist   = list(train = dtrain.f1, test = dval.f1),
      nrounds     = 20,
      eval_metric = "logloss",
      alpha       = 15,
      lambda      = 0,
      min_child_weight = 1,
      objective   = "binary:logistic"
    )
    
    loss.list[[omi]] <- ggplot(
      melt(model.f1$evaluation_log, id.vars = "iter"),
      aes(x = iter, y = value, color = variable)
    ) + geom_line() + theme_bw()
    
    # =========================================================================
    # STEP 4 – Run all 25 folds
    # =========================================================================
    
    # Pre-allocate score data frame (one row per fold)
    score_list <- data.frame(
      folds     = seq_along(names(folds)),
      scores    = 0,
      AUC       = 0,
      AUC.top10 = 0,
      rounds    = 0
    )
    
    for (i in names(folds)) {
      
      cat("Fold:", i, "\n")
      
      x                 <- folds[[i]]
      training_fold     <- .fix_colnames(train_data[-x, , drop = FALSE], omi, standards.variables)
      val_fold          <- .fix_colnames(train_data[x,  , drop = FALSE], omi, standards.variables)
      train_labels.fold <- train_labels[-x]
      val_labels.fold   <- train_labels[x]
      
      all_ids                <- row.names(train_data)
      train.labels.list[[i]] <- all_ids[-x]
      val.labels.list[[i]]   <- all_ids[x]
      
      training_xgb <- xgb.DMatrix(data = training_fold, label = train_labels.fold)
      testing_xgb  <- xgb.DMatrix(data = val_fold,      label = val_labels.fold)
      
      # ----- XGBoost model (full feature set) --------------------------------
      set.seed(seed)
      best_out <- xgboost(
        data                  = training_xgb,
        nround                = xgb_params$nround,
        alpha                 = alpha.omic,
        lambda                = xgb_params$lambda,
        eta                   = xgb_params$eta,
        gamma                 = xgb_params$gamma,
        max_depth             = xgb_params$max_depth,
        early_stopping_rounds = xgb_params$early_stopping_rounds,
        eval_metric           = xgb_params$eval_metric,
        objective             = xgb_params$objective,
        nthread               = xgb_params$nthread
      )
      
      best_rounds <- best_out$best_iteration
      fold_idx    <- which(i == names(folds))
      
      pred <- predict(best_out, testing_xgb,
                      iteration_range = c(1, best_rounds))
      
      auc.fold  <- Metrics::auc(actual = val_labels.fold, predicted = pred)
      pred.class <- ifelse(pred >= 0.5, 1, 0)
      
      model_list[[i]]            <- best_out
      labelsValidation.list[[i]] <- val_labels.fold
      predictions.list[[i]]      <- pred
      confusionMatrix.list[[i]]  <- confusionMatrix(
        factor(pred.class), factor(val_labels.fold)
      )
      score_list[fold_idx, "rounds"] <- best_rounds
      score_list[fold_idx, "scores"] <- best_out$evaluation_log[[2]][best_rounds]
      score_list[fold_idx, "AUC"]    <- auc.fold
      
      if (!omi %in% standards.variables) {
        importance_matrix            <- xgb.importance(model = best_out)
        feature.importance.list[[i]] <- importance_matrix
      }
      
      # =======================================================================
      # STEP 4.1 – Top-N feature sub-models (top 2 … top)
      # =======================================================================
      if (!omi %in% standards.variables) {
        
        for (u in 2:top) {
          
          feat_label <- paste0("top.", u)
          
          if (nrow(importance_matrix) >= u) {
            
            top_feats <- importance_matrix[
              order(importance_matrix$Gain, decreasing = TRUE), ]$Feature[1:u]
            
            dtrain.fs <- xgb.DMatrix(
              data  = training_fold[, top_feats, drop = FALSE],
              label = train_labels.fold
            )
            dval.fs <- xgb.DMatrix(
              data  = val_fold[, top_feats, drop = FALSE],
              label = val_labels.fold
            )
            
            set.seed(seed)
            model.fs <- xgb.train(
              data                  = dtrain.fs,
              watchlist             = list(train = dtrain.fs, eval = dval.fs),
              nrounds               = xgb_params$nround,
              alpha                 = alpha.omic,
              lambda                = xgb_params$lambda,
              eta                   = xgb_params$eta,
              gamma                 = xgb_params$gamma,
              max_depth             = xgb_params$max_depth,
              early_stopping_rounds = 10,
              eval_metric           = xgb_params$eval_metric,
              objective             = xgb_params$objective,
              nthread               = xgb_params$nthread,
              verbose               = 0
            )
            
            best_rounds.fs <- model.fs$best_iteration
            
            pred.fs <- predict(model.fs, dval.fs,
                               iteration_range = c(1, best_rounds.fs))
            
            auc.fs        <- Metrics::auc(actual = val_labels.fold,
                                          predicted = pred.fs)
            pred.class.fs <- ifelse(pred.fs >= 0.5, 1, 0)
            
            model_list.top[[feat_label]][[i]]              <- model.fs
            predictions.list.top[[feat_label]][[i]]        <- pred.fs
            confusionMatrix.list.top[[feat_label]][[i]]    <- confusionMatrix(
              factor(pred.class.fs), factor(val_labels.fold)
            )
            feature.importance.list.top[[feat_label]][[i]] <- xgb.importance(
              model = model.fs
            )
        
            #      store per-feature-count AUC in a dedicated named column
            col_name <- paste0("AUC.", feat_label)
            if (!col_name %in% colnames(score_list)) {
              score_list[[col_name]] <- NA_real_
            }
            score_list[fold_idx, col_name] <- auc.fs
            
          } else {
            model_list.top[[feat_label]][[i]]              <- NULL
            predictions.list.top[[feat_label]][[i]]        <- NULL
            confusionMatrix.list.top[[feat_label]][[i]]    <- NULL
            feature.importance.list.top[[feat_label]][[i]] <- NULL
            
            col_name <- paste0("AUC.", feat_label)
            if (!col_name %in% colnames(score_list)) {
              score_list[[col_name]] <- NA_real_
            }
            score_list[fold_idx, col_name] <- auc.fold
          }
        } # end top-N loop
      }   # end standards.variables guard
      
    } # end fold loop
    
    # =========================================================================
    # STEP 5 – Aggregate per-omic results into output lists
    # =========================================================================
    score_list.all[[omi]]              <- score_list
    predictions.list.all[[omi]]        <- predictions.list
    predictions.list.all.top[[omi]]    <- predictions.list.top
    model_list.all[[omi]]              <- model_list
    model_list.all.top[[omi]]          <- model_list.top
    feature.importance.list.all[[omi]] <- feature.importance.list
    feature.importance.list.top.all[[omi]] <- feature.importance.list.top
    confusionMatrix.list.all[[omi]]    <- confusionMatrix.list
    confusionMatrix.lit.top.all[[omi]] <- confusionMatrix.list.top
    train.labels.list.out[[omi]]       <- train.labels.list
    val.labels.list.out[[omi]]         <- val.labels.list
    
    # =========================================================================
    # STEP 5a – Pooled cross-validated AUC (full feature model)
    # =========================================================================
    cat("STEP 5: Computing pooled cvAUC for", omi, "\n")
    
    out    <- cvAUC(predictions.list, labelsValidation.list,
                    label.ordering = NULL)
    out.ci <- ci.pooled.cvAUC(predictions.list, labelsValidation.list,
                              label.ordering = NULL, folds = NULL,
                              confidence = 0.95,    ids = folds)
    
    
    
    out.list[[omi]]    <- out
    cvAUC.list[[omi]]  <- list(cvAUC = out, out.ci = out.ci)
    cvAUC.value[[omi]] <- out$cvAUC
    cvAUC.se[[omi]]    <- out.ci$se
    
    # =========================================================================
    # STEP 6 – Pooled cvAUC for each top-N sub-model
    # =========================================================================
    if (!omi %in% c(standards.variables, "risk", "Confounders")) {
      
      cat("STEP 6: Computing pooled cvAUC for top-N models –", omi, "\n")
      
      for (u in names(predictions.list.top)) {
        
        preds.u <- predictions.list.top[[u]]
        if (is.null(preds.u)) next
        
        preds.u  <- preds.u[!is.na(names(preds.u))]
        labels.u <- labelsValidation.list[names(preds.u)]
        
        for (l in names(preds.u)) {
          n_pred  <- length(preds.u[[l]])
          n_label <- length(labels.u[[l]])
          if (n_pred < n_label) {
            labels.u[[l]] <- labels.u[[l]][seq_len(n_pred)]
          } else if (n_pred > n_label) {
            preds.u[[l]]  <- preds.u[[l]][seq_len(n_label)]
          }
        }
        
        out.fs <- tryCatch(
          cvAUC(preds.u, labels.u, label.ordering = NULL),
          error = function(e) NA
        )
    
        out.ci.fs <- tryCatch(
          ci.pooled.cvAUC(preds.u, labels.u, label.ordering = NULL,
                          folds = NULL, confidence = 0.95,ids = folds),
          error = function(e) NA
        )
        
        out.list.top[[omi]][[u]]   <- out.fs
        cvAUC.list.top[[omi]][[u]] <- list(cvAUC = out.fs, out.ci = out.ci.fs)
      }
    }
    
  } # end omic loop
  
  # ---------------------------------------------------------------------------
  # Return everything
  # ---------------------------------------------------------------------------
  list(
    loss.list.all               = loss.list,
    cvAUC.list.all              = cvAUC.list,
    cvAUC.list.top.all          = cvAUC.list.top,
    feature.importance.list     = feature.importance.list.all,
    feature.importance.list.top = feature.importance.list.top.all,
    feature.importance.list2    = feature.importance.list2,
    cvAUC.value                 = cvAUC.value,
    cvAUC.se                    = cvAUC.se,
    model_list.all              = model_list.all,
    model_list.all.top          = model_list.all.top,
    predictions.list            = predictions.list.all,
    predictions.list.top        = predictions.list.all.top,
    labelsValidation.list       = labelsValidation.list,
    folds.all                   = folds.list,
    out                         = out.list,
    out.top                     = out.list.top,
    score_list.all              = score_list.all,
    confusionMatrix.list.all    = confusionMatrix.list.all,
    confusionMatrix.lit.top.all = confusionMatrix.lit.top.all,
    train.labels.ID.out         = train.labels.list.out,
    val.labels.ID.out           = val.labels.list.out
  )
}


# =============================================================================
# Internal helper – fix column names for standards variables
# (single-column matrices lose their name after sub-setting)
# =============================================================================
.fix_colnames <- function(mat, omi, standards.variables) {
  if (omi %in% standards.variables) {
    mat <- as.matrix(mat)
    colnames(mat)[1] <- omi
  }
  mat
}

