# Author: Camila Alvarez-Silva

library(dplyr)

source("/projects/arumugam/people/xhf865/Camila/GALAXY/GALA/xgboost/R_manuscript/Final/Cleaned/Github/xgboost_eval.GALAXY.R")



# Evaluate models with top 15 features
F1.bestmodels <- xgboost_eval(xgboost.object = xgboost.Fibrosis1, top = 15)
F2.bestmodels <- xgboost_eval(xgboost.object = xgboost.Fibrosis2, top = 15)
Inflammation.bestmodels <- xgboost_eval(xgboost.object = xgboost.Inflammation, top = 15)
steatosis.bestmodels <- xgboost_eval(xgboost.object = xgboost.steatosis, top = 15)

# Extract AUC results and assign panel names
F1.bestmodels2 <- F1.bestmodels$AUC.all.top
F1.bestmodels2$panel <- "F>=2"

F2.bestmodels2 <- F2.bestmodels$AUC.all.top
F2.bestmodels2$panel <- "F>=3"

Inflammation.bestmodels2 <- Inflammation.bestmodels$AUC.all.top
Inflammation.bestmodels2$panel <- "I>=2"

steatosis.bestmodels2 <- steatosis.bestmodels$AUC.all.top
steatosis.bestmodels2$panel <- "S>=2"

# Combine all models into a single dataset
bestmodels <- rbind(F1.bestmodels2, F2.bestmodels2, Inflammation.bestmodels2, steatosis.bestmodels2)

# Remove rows with missing AUC values
bestmodels <- bestmodels %>% filter(!is.na(AUC))

bestmodels$AUC <- round(bestmodels$AUC, 2)
bestmodels$sd <- round(bestmodels$sd, 2)

# Rename columns for clarity
colnames(bestmodels)[2] <- "AUC.cv"
colnames(bestmodels)[4] <- "features"
bestmodels$no.features <- as.numeric(str_remove_all(bestmodels$features, "top."))

# Group by 'panel' and 'omic', then find the maximum AUC.cv and corresponding features
bestmodels.auc <- bestmodels %>% 
  group_by(panel, omic) %>% 
  summarise(
    max.AUC.cv = max(AUC.cv),
    features = features[which(AUC.cv == max(AUC.cv))] 
  )

# Extract numeric feature count from the 'features' column 
bestmodels.auc$no.features <- as.numeric(str_remove_all(bestmodels.auc$features, "top."))

# Since multiple  models can have the same maximun AUC we filter to retain only models with the minimum number of features within each 'panel' and 'omic' group
bestmodels.auc.filt <- bestmodels.auc %>% 
  group_by(panel, omic) %>% 
  filter(no.features == min(no.features))

# Compute 99% of the max AUC.cv and round it to two decimal places
bestmodels.auc.filt$AV.99 <- round(bestmodels.auc.filt$max.AUC.cv * 0.99, 2)

# Initialize an empty list to store the best AUC models
best.auc.list <- list()

# Iterate over each row in the filtered bestmodels.auc.filt dataframe
for (i in 1:nrow(bestmodels.auc.filt)) {
  df <- bestmodels.auc.filt[i, ]
  panel <- df$panel
  omic <- as.character(df$omic)
  
  # Filter original data for the current panel and omic
  df.filt <- bestmodels %>% 
    filter(panel == df$panel & omic == df$omic)
  
  #  Filter rows with AUC.cv greater than or equal to 99% of max AUC.cv
  df.filt2 <- df.filt %>%  filter(AUC.cv >= df$AV.99)
  
  # Keep only the rows with the minimum number of features
  df.filt2 <- df.filt2 %>%  filter(no.features == min(no.features))
  
  # Store the filtered data in a nested list structure
  best.auc.list[[panel]][[omic]] <- df.filt2
}

# Convert nested list to data frames for each panel
best.auc.df.f2 <- do.call("rbind", best.auc.list$`F>=2`)
best.auc.df.f3 <- do.call("rbind", best.auc.list$`F>=3`)
best.auc.df.I <- do.call("rbind", best.auc.list$`I>=2`)
best.auc.df.S <- do.call("rbind", best.auc.list$`S>=2`)

# Combine results into a final dataframe
best.auc.df <- rbind(best.auc.df.f2, best.auc.df.f3, best.auc.df.I, best.auc.df.S)

saveRDS("best.auc.df")