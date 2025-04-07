# Author: Camila Alvarez-Silva

source("/../../xgboost_eval.GALAXY.R")

# Load test dataset
df.matrix.test<-readRDS("/../../df.matrix.test.rds")

# Load test metadata

metadata<-as.data.frame(fread("/../../Camila_test_data.2025.csv"))
row.names(metadata)<-metadata$SampleID

# Load best models (Result from xgboost best model selection)
bestmodels.filt.AUC<-readRDS("/../..//bestmodels.filt.AUC")

# Load trained XGBoost model (Result from xgboost_train)
xgboost.Fibrosis1<-readRDS("/../xgboost.Fibrosis1")


# Select features of interest from test metadata
features <- metadata %>% dplyr::select("fibrosis", "kleiner_numeric", "inflam_numeric", "steatosis_numeric")

# Transform feature values into binary classifications
features$fibrosis <- ifelse(features$kleiner_numeric >= 2, 1, 0)
features$kleiner_numeric <- ifelse(features$kleiner_numeric >= 3, 1, 0)
features$inflam_numeric <- ifelse(as.numeric(features$inflam_numeric) < 2, 0, 1)
features$steatosis_numeric <- ifelse(as.numeric(features$steatosis_numeric) < 2, 0, 1)

colnames(features) <- c("F>=2", "F>=3", "I>=2", "S>=2")

# Function to validate the XGBoost model on test data


test.validation=function(xgboost.result=xgboost.result, panel1=panel, df.test=df.test, features=features, top.features=top.features,bestmodels.filt.AUC=bestmodels.filt.AUC, metadata.test=metadata){
  
  
  #	Identifies the best model per omic reduced model across all folds (25)
  
  bestmodels<-xgboost_eval(xgboost.result, top=15)
  bestmodels.df<-bestmodels$bestModel
  bestmodels.df$omic.feature<-paste(bestmodels.df$omic,bestmodels.df$features, sep = ".")
  
  #bestmodels$AUC.barplot.all
  
  bestmodels.filt.AUC$panel<-as.character(bestmodels.filt.AUC$panel)
  bestmodels.df.top<-bestmodels.filt.AUC %>% dplyr::filter(panel==panel1)
  bestmodels.df.top$omic.feature<-paste(bestmodels.df.top$omic,bestmodels.df.top$features, sep = ".")
  
  bestmodels.filt<-bestmodels.df %>% filter(omic.feature %in% bestmodels.df.top$omic.feature)
  
  
  df.cofounders.risk<-data.frame("omic"= c("Confounders","risk"), 
                                 "features"= c(NA,NA), 
                                 "maxAUC"= c(NA,NA), 
                                 "model"=c(NA,NA), 
                                 "fold.name"= c(NA,NA), 
                                 "AUC.cv"= c(NA,NA), 
                                 "sd"= c(NA,NA), 
                                 "no.features" =c(NA,NA), 
                                 "omic.feature"=c(NA,NA))
  
  
  #Since we did not do feature reducction for eather cofounder or genetic resiscore  have to  add the information here
  ### Confounders
  df.cofounders.risk$features[[1]]<-"top.11"
  df.cofounders.risk$maxAUC[[1]]<-max(xgboost.result$cvAUC.list$Confounders$cvAUC$fold.AUC)
  max.auc<-which((xgboost.result$cvAUC.list$Confounders$cvAUC$fold.AUC)==df.cofounders.risk$maxAUC[[1]])
  df.cofounders.risk$model[[1]]<-max.auc
  df.cofounders.risk$fold.name[[1]]<-names(xgboost.result$model_list.all$Confounders)[max.auc]
  df.cofounders.risk$AUC.cv[[1]]<-xgboost.result$cvAUC.list$Confounders$cvAUC$cvAUC
  df.cofounders.risk$sd[[1]]<-sd(xgboost.result$cvAUC.list$Confounders$cvAUC$fold.AUC)
  df.cofounders.risk$no.features[[1]]<-11
  df.cofounders.risk$omic.feature[[1]]<-"Confounders.top.11"
  
  ### risk
  df.cofounders.risk$features[[2]]<-"top.2"
  df.cofounders.risk$maxAUC[[2]]<-max(xgboost.result$cvAUC.list$risk$cvAUC$fold.AUC)
  max.auc<-which((xgboost.result$cvAUC.list$risk$cvAUC$fold.AUC)==df.cofounders.risk$maxAUC[[2]])
  df.cofounders.risk$model[[2]]<-max.auc
  df.cofounders.risk$fold.name[[2]]<-names(xgboost.result$model_list.all$risk)[max.auc]
  df.cofounders.risk$AUC.cv[[2]]<-xgboost.result$cvAUC.list$risk$cvAUC$cvAUC
  df.cofounders.risk$sd[[2]]<-sd(xgboost.result$cvAUC.list$risk$cvAUC$fold.AUC)
  df.cofounders.risk$no.features[[2]]<-2
  df.cofounders.risk$omic.feature[[2]]<-"risk.top.2"
  
  bestmodels.filt<-rbind(bestmodels.filt,df.cofounders.risk)
  bestmodels.filt$panel<-panel1
  bestmodels.filt$test.AUC<-NA
  bestmodels.filt$dim.test<-NA
  bestmodels.filt$dim.Disease<-NA
  
  #extract features  of each selected model
  topFeatures<-reshape2::melt(bestmodels$topFeatures.list)
  colnames(topFeatures)<-c("Feature","n.Features", "omic")
  topFeatures$omic.feature<-paste(topFeatures$omic,topFeatures$n.Features, sep = ".")
  topFeatures$panel<-panel1
  topFeatures<-topFeatures %>% filter(omic.feature %in% bestmodels.df.top$omic.feature)
  
  
  #select target of interest to validate on the test metadata
  trait<-features %>% dplyr::select(as.name(panel1))
  trait<-trait %>% dplyr::filter(!is.na(trait[,1]))
  
  confusionMatrix.list<-list()
  df.prediction.list<-list()
  metadata.omic.test.list<-list()
  ####################################
  # Try model
  #################################### 
  
  for(i in 1:nrow(bestmodels.filt)){
    
    omic<-bestmodels.filt$omic[i]
    
    if(omic %in% c( "Confounders", "risk")){
      model.names<-names(xgboost.result$model_list.all[[omic]])
      model<-model.names[bestmodels.filt$model[i]]
      
      omic.feature.model<-bestmodels.filt$omic.feature[i]
      
      model.xgboost<-xgboost.result$model_list.all[[omic]][[model]]
      best_rounds<- which.min(model.xgboost$evaluation_log[[2]])}else{
        
        feature<-as.character(bestmodels.filt$features[i])
        model.names<-names(xgboost.result$model_list.all.top[[omic]][[feature]])
        model<-model.names[bestmodels.filt$model[i]]
        
        omic.feature.model<-bestmodels.filt$omic.feature[i]
        
        model.xgboost<-xgboost.result$model_list.all.top[[omic]][[feature]][[model]]
        best_rounds<- which.min(model.xgboost$evaluation_log[[2]])}
    
    #Filter df.matrix by selected features
    features.model<-model.xgboost$feature_names
    
    test<-df.test %>% dplyr::select(all_of(features.model))
    
    test <- test[rowSums(is.na(test)) < ncol(test), ]
    trait2<- trait[rownames(test),]
    test = as.matrix(test)
        test_label=trait2
    
    
    # put our test data into two separates Dmatrixs objects
    dtest <- xgb.DMatrix(data=test, label= test_label)
    
    #use the selected model to evalualte on the dtest
    pred.fs =predict(model.xgboost, dtest,iteration_range = best_rounds)
    auc.test<-Metrics::auc(actual=test_label, predicted= pred.fs)
    
    #confusionMatrix 
    pred.fact.fs<-ifelse(pred.fs >=0.5,1,0)
    
    # 
    #DF predictions
    df<-data.frame(ID=rownames(test),trait=test_label,Prediction=pred.fs,omic=omic,panel1=panel1, No.samples=nrow(test))
    
    
    #Add disease number of patients 
    metadata.omic.test<-filter(metadata,SampleID %in% row.names(test))
    features.omic<-features[row.names(test),panel1]
    disease<-table(features.omic)[2]
    bestmodels.filt$test.AUC[i]<-auc.test
    bestmodels.filt$dim.test[i]<-nrow(test)
    bestmodels.filt$dim.Disease[i]<-disease
    bestmodels.filt$fold.name[i]<-model
    
    
    df.prediction.list[[panel1]][[omic]]<-df
    metadata.omic.test.list[[panel1]][[omic]]<-metadata.omic.test
  }
  
  return.list<-list(bestmodels.filt=bestmodels.filt,
                    metadata.omic.test=metadata.omic.test.list,
                    df.prediction.list=df.prediction.list)
  return(return.list)
}





# Run validation for Fibrosis
Fibrosis1 <- test.validation(
  xgboost.result = xgboost.Fibrosis1,
  panel1 = "F>=2",
  df.test = df.matrix.test,
  features = features,
  bestmodels.filt.AUC = bestmodels.filt.AUC,
  metadata.test = metadata)
