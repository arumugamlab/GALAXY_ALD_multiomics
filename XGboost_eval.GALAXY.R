xgboost_eval <- function(xgboost.object, top) {
  
  # Function to generate performance plots for an XGBoost model
  # INPUT:
  #   - xgboost_object: A trained XGBoost model with AUC values and feature importance
  #   - top: Number of top features to consider
  # OUTPUT:
  #   - A list of plots and performance metrics
  
  
  # Initialize lists to store outputs
  plot.AUC.all <- list()
  cvAUC.list <- list()
  topFeatures.list <- list()
  upset.list <- list()
  plotAUC.list <- list()
  plot.complete <- list()
  confusionMatrix.list <- list()
  omics <- names(xgboost.object$cvAUC.value)
  cvAUC.df.top.list <- list()
  
  
  # Define color palette for plots
  color_palette <- c(
    "#F3DFC1", "#E76F51", "#F4A261", "#E9C46A", "#2A9D8F", "#457B9D", "#264653", "#A8DADC", "#CAD2C5", "#84A98C",
    "#52796F", "#354F52", "#6D6875", "#B5838D", "#E5989B", "#FFB4A2", "#FFCDB2", "#8E7DBE", "#A4A6D0", "#B3C5F1",
    "#D0E1FD", "#CDEFF7", "#7FD3C3", "#5EB5A6", "#3C848D", "#2F5E6E", "#4D4C75", "#5D5A89", "#6E7AA7", "#83A1C7",
    "#A7BCD6", "#C1D2E6", "#DDE5F3", "#E8EDF7", "#FACFD2", "#E0B1B8", "#C2979E", "#A67E84", "#8C6368", "#724C4F",
    "#5E4043", "#4C3A39", "#352F2F", "#2C2323", "#231919", "#8A8686", "#A7A5A5", "#C3C1C1", "#E0DFDF", "#F0F0F0", "#F3DFC1"
  )
  
  # --- Generate AUC Bar Plots: ---
  
  
  # Extract cross-validation AUC values
  cvAUC.value <- xgboost.object$cvAUC.value
  cvAUC.df <- reshape2::melt(cvAUC.value)
  cvAUC.df <- cvAUC.df[order(cvAUC.df$value, decreasing = TRUE), ]
  cvAUC.df$L1 <- factor(cvAUC.df$L1, levels = cvAUC.df$L1)
  
  # Compute standard deviation of AUC values
  cvAUC.list <- xgboost.object$cvAUC.list.all   #xgboost.object$cvAUC.list.all 
  sdAUC.val.df <- data.frame(L1 = names(cvAUC.list), sd = NA)
  
  
  for (i in seq_along(cvAUC.list)) {
    sdAUC.val.df$sd[i] <- sd(cvAUC.list[[i]]$cvAUC$fold.AUC)
  }
  
  # Merge standard deviation values
  cvAUC.df <- merge(cvAUC.df, sdAUC.val.df, by = "L1")
  cvAUC.df <- cvAUC.df[order(cvAUC.df$value, decreasing = FALSE), ]
  cvAUC.df$L1 <- factor(cvAUC.df$L1, levels = cvAUC.df$L1)
  
  
  # Generate AUC bar plot for all omics
  AUC.barplot.all.omics <- ggplot(data = cvAUC.df, aes(x = value, y = L1)) +
    geom_bar(stat = "identity", width = 0.6, alpha = 0.9, fill = "#264653") +
    theme_bw() +
    theme(legend.position = "none") +
    geom_text(aes(label = round(value, 2)), hjust = 1.1, color = "white") +
    xlab("average AUC") +
    ylab(NULL) +
    geom_errorbar(aes(xmin = (value - (sd)), xmax = (value + (sd))), width = .2, position = position_dodge(.9), color = "grey")
  
  # --- TOP FEATURES PLOTTING ---
  
  # Extracts performance metrics for models across top features models
  # Generates AUC plots with error bars and labels showing top feature performance.
  
  
  # Process and plot AUC across each  top feature models 
  
  top.df <- data.frame(omic = names(xgboost.object$cvAUC.list.top.all), top = NA)
  rownames(top.df)<-top.df$omic
  
  for (j in names(xgboost.object$cvAUC.list.top.all)) {
    top.n <- (length(xgboost.object$cvAUC.list.top.all[[j]]) + 1)
    
    top.df[j,2] <- top.n
  }
  
  top.max <- top.df$omic[which(top.df$top == top)]
  
  for (j in 1:(top - 1)) {
    cvAUC.df.top <- data.frame(omic = names(xgboost.object$cvAUC.list.top.all), AUC = NA, sd = NA)
    rownames(cvAUC.df.top)<-cvAUC.df.top$omic
    for (i in names(xgboost.object$cvAUC.list.top)) {
      omi <- i
      
      if (omi %in% top.max) {
        cvAUC.list <- xgboost.object$cvAUC.list.top.all[[i]][[j]]
        
        if (!anyNA(cvAUC.list)) {
          AUC.top <- cvAUC.list$cvAUC$cvAUC
          cvAUC.df.top[i, 2] <- AUC.top
          
          sdAUC.val <- sd(cvAUC.list$cvAUC$fold.AUC)
          if (!is.na(sdAUC.val)) {
            cvAUC.df.top[i, 3] <- sdAUC.val
          } else {
            cvAUC.df.top[i, 3] <- 0
          }
        }
        
        features <- names(xgboost.object$cvAUC.list.top[[i]][j])
        
        # Sort and plot AUC per omic
        cvAUC.df.top <- cvAUC.df.top[order(cvAUC.df.top$AUC, decreasing = FALSE), ]
        cvAUC.df.top$omic <- factor(cvAUC.df.top$omic, levels = cvAUC.df.top$omic)
        
        p <- ggplot(data = cvAUC.df.top, aes(x = AUC, y = omic, fill = omic)) +
          geom_bar(stat = "identity", width = 0.6, alpha = 0.9) +
          theme_bw() +
          scale_fill_manual(values = color_palette) +
          theme(legend.position = "none") +
          geom_text(aes(label = round(AUC, 2)), hjust = 1.1, color = "white") +
          xlab("average AUC [top features model]") +
          ylab(NULL) +
          ggtitle(features) +
          geom_errorbar(aes(xmin = (AUC - (sd)), xmax = (AUC + (sd))), width = .2, position = position_dodge(.9), color = "grey")
        
        plot.AUC.all[[features]] <- p
        cvAUC.df.top.list[[features]] <- cvAUC.df.top
      } else {
        omi <- names(xgboost.object$cvAUC.list.top)[i]
        top.n <- top.df[i,2]
        if (j <= (top.n - 1)) {
          cvAUC.list <- xgboost.object$cvAUC.list.top[[i]][[j]]
          
          if (!anyNA(cvAUC.list)) {
            AUC.top <- cvAUC.list$cvAUC$cvAUC
            cvAUC.df.top[i, 2] <- AUC.top
            
            sdAUC.val <- sd(cvAUC.list$cvAUC$fold.AUC)
            if (!is.na(sdAUC.val)) {
              cvAUC.df.top[i, 3] <- sdAUC.val
            } else {
              cvAUC.df.top[i, 3] <- 0
            }
          }
          
          features <- names(xgboost.object$cvAUC.list.top[[i]][j])
          
          cvAUC.df.top <- cvAUC.df.top[order(cvAUC.df.top$AUC, decreasing = FALSE), ]
          cvAUC.df.top$omic <- factor(cvAUC.df.top$omic, levels = cvAUC.df.top$omic)
          
          p <- ggplot(data = cvAUC.df.top, aes(x = AUC, y = omic, fill = omic)) +
            geom_bar(stat = "identity", width = 0.6, alpha = 0.9) +
            theme_bw() +
            scale_fill_manual(values = color_palette) +
            theme(legend.position = "none") +
            geom_text(aes(label = round(AUC, 2)), hjust = 1.1, color = "white") +
            xlab("average AUC [top features model]") +
            ylab(NULL) +
            ggtitle(features) +
            geom_errorbar(aes(xmin = (AUC - (sd)), xmax = (AUC + (sd))), width = .2, position = position_dodge(.9), color = "grey")
          
          plot.AUC.all[[features]] <- p
          cvAUC.df.top.list[[features]] <- cvAUC.df.top
        }
      }
    }
  }
  
  
  
  # Generate table with performance of feature reduced models across the X folds
  
  cvAUC.list <- list()
  
  for (j in names(xgboost.object$cvAUC.list.top.all)) {
    top <- (length(xgboost.object$cvAUC.list.top.all[[j]]) + 1)
    
    for (i in 1:(top - 1)) {
      cvAUC.df <- xgboost.object$cvAUC.list.top.all[[j]][[i]]
      if (!anyNA(cvAUC.df)) {
        n.folds <- length(cvAUC.df$cvAUC$fold.AUC)
        
        cvAUC.df.top <- data.frame(omic = rep(NA, n.folds), features = rep(NA, n.folds), AUC = rep(NA, n.folds), AUC.cv = rep(NA, n.folds), sd = rep(NA, n.folds))
        feature <- names(xgboost.object$cvAUC.list.top.all[[j]][i])
        omic <- names(xgboost.object$cvAUC.list.top.all[j])
        
        cvAUC.df <- xgboost.object$cvAUC.list.top.all[[j]][[i]]
        
        AUC.top <- cvAUC.df$cvAUC$fold.AUC
        sdAUC.val <- sd(AUC.top)
        AUC.cv <- cvAUC.df$cvAUC$cvAUC
        cvAUC.df.top$omic <- omic
        cvAUC.df.top$features <- feature
        cvAUC.df.top$AUC <- AUC.top
        cvAUC.df.top$AUC.cv <- AUC.cv
        cvAUC.df.top$sd <- sdAUC.val
        to.return <- list()
        cvAUC.list[[omic]][[feature]] <- cvAUC.df.top
      }
    }
  }
  
  cvAUC.list.melt <- reshape2::melt(cvAUC.list, id.vars = c("omic", "features", "AUC", "AUC.cv", "sd"))
  
  cvAUC.list.melt$L2 <- NULL
  cvAUC.list.melt$L1 <- NULL
  
  top.max <- max(as.numeric(str_split_fixed(unique(cvAUC.list.melt$features), "\\.", 2)[, 2]))
  cvAUC.list.melt$features <- factor(cvAUC.list.melt$features, levels = paste0("top.", 2:top.max))
  
  cvAUC.list.melt2 <- cvAUC.list.melt
  cvAUC.list.melt2$AUC <- NULL
  cvAUC.list.melt2 <- unique(cvAUC.list.melt2)
  cvAUC.list.melt2$AUC.cv <- round(cvAUC.list.melt2$AUC.cv, 2)
  cvAUC.list.melt2$features <- str_remove_all(cvAUC.list.melt2$features, "top.")
  cvAUC.list.melt2$features <- factor(cvAUC.list.melt2$features, levels = c(2:top.max))
  cvAUC.list.melt2$sd[which(is.na(cvAUC.list.melt2$sd) == T)] <- 0
  
  # --- AUC curves across all the feature reduced models for comparison ---
  
  AUC.top.feature <- ggplot(cvAUC.list.melt2, aes(x = features, y = AUC.cv, group = 1, color = omic, label = AUC.cv)) +
    geom_line() +
    theme_bw() +
    theme(legend.position = "none") +
    facet_wrap(~omic, ncol = 4, scales = "free") +
    scale_color_manual(values = color_palette) +
    geom_point() +
    ylab("cv.AUC") +
    geom_label_repel(box.padding = 0.35, point.padding = 0.5, size = 2) +
    geom_errorbar(aes(ymin = (AUC.cv - sd), ymax = (AUC.cv + sd)), width = .2, position = position_dodge(0.05))
  
  
  # --- Best Model Identification
  
  #	Identifies the best model per omic reduced model based on maximum AUC an minimun number of features.
  #	Extracts feature importance tables for each best-performing model.
  
  
  # --- Extract top-performing model information  ---
  
  
  bestModel <- cvAUC.list.melt %>%
    group_by(omic, features) %>%
    dplyr::summarise(
      maxAUC = max(AUC),
      model = which(AUC == max(AUC)),
      AUC.cv = unique(AUC.cv),
      sd = unique(sd)
    )
  
  for (i in 1:nrow(bestModel)) {
    omic <- bestModel$omic[i]
    feature <- as.character(bestModel$features[i])
    model <- bestModel$model[i]
    maxAUC <- bestModel$maxAUC[i]
    df.features <- as.data.frame(xgboost.object$feature.importance.list.top[[omic]][[feature]][model])
    colnames(df.features) <- c("Feature", "Gain", "Cover", "Frequency")
    
    topFeatures <- df.features$Feature
    
    topFeatures.list[[omic]][[feature]] <- topFeatures
    
    # --- Confusion Matrix Visualization---
    
    # Confusion matrix for the best-performing model.
    
    confusionMatrix.bm <- data.frame(xgboost.object$confusionMatrix.lit.top.all[[omic]][[feature]][[model]]$table)
    
    plotTable <- confusionMatrix.bm %>%
      mutate(goodbad = ifelse(confusionMatrix.bm$Prediction == confusionMatrix.bm$Reference, "good", "bad")) %>%
      group_by(Reference) %>%
      mutate(prop = Freq / sum(Freq))
    
    confusionMatrix <- ggplot(data = plotTable, mapping = aes(x = Reference, y = Prediction, fill = goodbad, alpha = prop)) +
      geom_tile() +
      geom_text(aes(label = Freq), vjust = .5, fontface = "bold", alpha = 1) +
      scale_fill_manual(values = c(good = "#264653", bad = "#2a9d8f")) +
      theme_bw() +
      theme(legend.position = "none") +
      xlim(rev(levels(confusionMatrix.bm$Reference))) +
      ggtitle(paste("Model:", feature, "Fold:", model, sep = " "))
    
    confusionMatrix.list[[omic]][[feature]] <- confusionMatrix
  }
  
  # --- Generate best model plot: higher AUC best lower number of features ---
  
  bestModel$no.features <- as.numeric(sapply(str_split(bestModel$features, "\\."), tail, 1))
  bestModel$AUC.cv <- round(bestModel$AUC.cv, 2)
  bestModel.top <- bestModel %>%
    group_by(omic) %>%
    dplyr::summarise(best = bestModel$features[which(AUC.cv == max(AUC.cv))])
  
  bestModel.top$no.features <- as.numeric(sapply(str_split(bestModel.top$best, "\\."), tail, 1))
  
  bestModel.top2 <- bestModel.top %>%
    group_by(omic) %>%
    dplyr::summarise(min = min(no.features))
  bestModel.top2$min <- paste0("top.", bestModel.top2$min)
  bestModel.top2$omic.feature <- paste(bestModel.top2$omic, bestModel.top2$min, sep = ".")
  
  AUC.all.melt <- reshape2::melt(cvAUC.df.top.list, id.vars = c("omic", "AUC", "sd"))
  AUC.all.melt$omic.feature <- paste(AUC.all.melt$omic, AUC.all.melt$L1, sep = ".")
  AUC.all.melt.filt <- AUC.all.melt %>% filter(AUC.all.melt$omic.feature %in% bestModel.top2$omic.feature)
  
  AUC.all.melt.filt <- AUC.all.melt.filt[order(AUC.all.melt.filt$AUC, decreasing = FALSE), ]
  AUC.all.melt.filt$omic.feature <- factor(AUC.all.melt.filt$omic.feature, levels = AUC.all.melt.filt$omic.feature)
  
  AUC.barplot.omics.all.bestmodels <- ggplot(data = AUC.all.melt.filt, aes(x = AUC, y = omic.feature, fill = omic)) +
    geom_bar(stat = "identity", width = 0.6, alpha = 0.9) +
    theme_bw() +
    scale_fill_manual(values = color_palette) +
    theme(legend.position = "none") +
    geom_text(aes(label = round(AUC, 2)), hjust = 1.1, color = "white") +
    xlab("average AUC [top features model]") +
    ylab(NULL) +
    geom_errorbar(aes(xmin = (AUC - (sd)), xmax = (AUC + (sd))), width = .2, position = position_dodge(.9), color = "grey")
  
  
  # --- Upset Plot Generation---
  
  # --- Plot cvAUC across feature reduced models including features included in each model---
  
  topFeatures.list.melt <- reshape2::melt(topFeatures.list)
  colnames(topFeatures.list.melt) <- c("Feature", "nFeatures", "omic")
  
  topFeatures.list.melt$Feature <- str_remove(topFeatures.list.melt$Feature, "_Relative Amount")
  topFeatures.list.melt$Feature <- str_remove(topFeatures.list.melt$Feature, " species incertae sedis")
  topFeatures.list.melt$Feature <- str_remove(topFeatures.list.melt$Feature, "mOTU_v25_")
  
  omics <- unique(topFeatures.list.melt$omic)
  
  
  for (i in omics) {
    test <- topFeatures.list.melt %>% filter(omic == i)
    test <- as.data.frame(t(unclass(table(test$Feature, test$nFeatures))))
    
    upsetPlot <- ComplexUpset::upset(test, intersect = colnames(test), base_annotations = list(), sort_intersections = "ascending", sort_intersections_by = "degree", set_sizes = FALSE)
    upset.list[[i]] <- upsetPlot
    
    cvAUC.list.melt.omic <- cvAUC.list.melt2 %>% filter(omic == i)
    
    p2 <- ggplot(cvAUC.list.melt.omic, aes(x = features, y = AUC.cv, group = 1, color = omic, label = AUC.cv)) +
      geom_line() +
      theme_bw() +
      theme(legend.position = "none") +
      scale_color_manual(values = color_palette) +
      geom_point() +
      ylab("cv.AUC") +
      geom_label_repel(box.padding = 0.35, point.padding = 0.5, size = 3) +
      ggtitle(i)
    plotAUC.list[[i]] <- p2
    
    p <- (p2 / upsetPlot)
    plot.complete[[i]] <- p
  }
  
  
  # ---------------------------------------------------------------------------
  # --- Optimal Feature-Reduced Model Selection (>= 99% of max AUC) ----------
  # ---------------------------------------------------------------------------
  #
  # Starting from AUC.all.melt (AUC.all.top), identifies the most parsimonious
  # model per panel x omic combination whose average AUC is >= 99% of the
  # maximum average AUC observed for that group.
  #
  # Inputs (from upstream steps):
  #   AUC.all.melt : melted data frame produced from cvAUC.df.top.list;
  #                  columns include omic, AUC, sd, L1 (feature level), omic.feature
  #
  # Outputs:
  #   best.auc.list : nested list [[panel]][[omic]] of the optimal model row(s)
  #   best.auc.df   : final data frame – one row per panel x omic, columns:
  #                     omic, AUC.cv, sd, panel, no.features
  # ---------------------------------------------------------------------------
  
  # Work from AUC.all.melt (exported as AUC.all.top); rename for clarity
  bestmodels <- AUC.all.melt
  
  # Remove rows with missing AUC values
  bestmodels <- bestmodels %>% filter(!is.na(AUC))
  bestmodels$AUC <- round(bestmodels$AUC, 2)
  bestmodels$sd  <- round(bestmodels$sd,  2)
  
  # Rename columns for clarity
  colnames(bestmodels)[colnames(bestmodels) == "AUC"] <- "AUC.cv"
  colnames(bestmodels)[colnames(bestmodels) == "L1"]  <- "features"
  
  # Derive numeric feature count from "top.X" labels
  bestmodels$no.features <- as.numeric(str_remove_all(bestmodels$features, "top\\."))
  
  # NOTE: AUC.all.melt does not carry a 'panel' column at this stage.
  # A 'panel' column is added here as a placeholder so the nested list
  # structure [[panel]][[omic]] is preserved. When xgboost_eval() is called
  # per clinical panel (steatosis, inflammation, fibrosis >=2, fibrosis >=3)
  # the caller should set xgboost.object$panel before passing in, or the
  # panel label can be injected by the wrapper script (xgboost_best_model).
  if (!is.null(xgboost.object$panel)) {
    bestmodels$panel <- xgboost.object$panel
  } else {
    bestmodels$panel <- "unknown"
  }
  
  # --- Step 1: find max AUC.cv per panel x omic and candidate features -------
  bestmodels.auc <- bestmodels %>%
    group_by(panel, omic) %>%
    dplyr::summarise(
      max.AUC.cv = max(AUC.cv),
      features   = features[which(AUC.cv == max(AUC.cv))],
      .groups    = "drop"
    )
  
  # Numeric feature count from "top.X" labels
  bestmodels.auc$no.features <- as.numeric(str_remove_all(bestmodels.auc$features, "top\\."))
  
  # --- Step 2: keep only the model with the fewest features at max AUC -------
  bestmodels.auc.filt <- bestmodels.auc %>%
    group_by(panel, omic) %>%
    dplyr::filter(no.features == min(no.features)) %>%
    dplyr::ungroup()
  
  # --- Step 3: compute 99% threshold per panel x omic row -------------------
  bestmodels.auc.filt$AV.99 <- round(bestmodels.auc.filt$max.AUC.cv * 0.99, 2)
  
  # --- Step 4: build best.auc.list -------------------------------------------
  best.auc.list <- list()
  
  for (i in 1:nrow(bestmodels.auc.filt)) {
    
    df    <- bestmodels.auc.filt[i, ]
    panel <- df$panel
    omic  <- as.character(df$omic)
    
    # Filter original data for the current panel and omic
    df.filt <- bestmodels %>%
      dplyr::filter(panel == df$panel & omic == df$omic)
    
    # Keep only rows with AUC.cv >= 99% of the maximum AUC.cv
    df.filt2 <- df.filt %>%
      dplyr::filter(AUC.cv >= df$AV.99)
    
    # Among those, keep only the rows with the minimum number of features
    df.filt2 <- df.filt2 %>%
      dplyr::filter(no.features == min(no.features))
    
    # Store in nested list [[panel]][[omic]]
    best.auc.list[[panel]][[omic]] <- df.filt2
  }
  
  # --- Step 5: flatten to best.auc.df ----------------------------------------
  best.auc.df <- do.call(rbind, lapply(names(best.auc.list), function(p) {
    do.call(rbind, lapply(names(best.auc.list[[p]]), function(o) {
      row <- best.auc.list[[p]][[o]]
      data.frame(
        omic        = o,
        AUC.cv      = row$AUC.cv,
        sd          = row$sd,
        panel       = p,
        no.features = row$no.features,
        stringsAsFactors = FALSE
      )
    }))
  }))
  
  rownames(best.auc.df) <- NULL
  
  
  # ---------------------------------------------------------------------------
  # Return results
  # ---------------------------------------------------------------------------
  
  to.return <- list(
    plot.complete                  = plot.complete,
    AUC.top.feature.curve          = AUC.top.feature,
    plot.AUC.all.feature.reduced   = plot.AUC.all,
    AUC.barplot.all.omics          = AUC.barplot.all.omics,
    confusionMatrix.best.model     = confusionMatrix.list,
    bestModel                      = bestModel,
    topFeatures.list               = topFeatures.list,
    AUC.barplot.omics.all.bestmodels = AUC.barplot.omics.all.bestmodels,
    bestModel.top                  = bestModel.top2,
    AUC.all.top                    = AUC.all.melt,
    cvAUC.list.melt                = cvAUC.list.melt2,
    best.auc.list                  = best.auc.list,
    best.auc.df                    = best.auc.df
  )
  return(to.return)
}
