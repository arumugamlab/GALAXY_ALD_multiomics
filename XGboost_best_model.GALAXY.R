# Author: Camila Alvarez-Silva

library(dplyr)

source("/../xgboost_eval.GALAXY.R")

#load xgboost models (output of XGboost_train)
xgboost.Fibrosis1<-readRDS("/../../xgboost.Fibrosis1.rds")
xgboost.Fibrosis2<-readRDS("/../../xgboost.Fibrosis2.rds")
xgboost.Inflammation<-readRDS("/../../xgboost.Inflammation.rds")
xgboost.Inflammation<-readRDS("/../../xgboost.Inflammation.rds")


# Evaluate models
xgboost.Fibrosis1$panel <- "F>=2"
F1.bestmodels <- xgboost_eval(xgboost.object = xgboost.Fibrosis1, top = 15)

xgboost.Fibrosis2$panel <- "F>=3"
F.bestmodels <- xgboost_eval(xgboost.object = xgboost.Fibrosis2, top = 15)

xgboost.Inflammation$panel <- "I>=2"
Inflammation.bestmodels <- xgboost_eval(xgboost.object = xgboost.Inflammation, top = 15)

xgboost.steatosis$panel <- "I>=2"
steatosis.bestmodels <- xgboost_eval(xgboost.object = xgboost.steatosis, top = 15)

# Combine all models into a single list
best.auc.list <- list(F1.bestmodels$best.auc.list, 
                      F.bestmodels$best.auc.list, 
                      Inflammation.bestmodels$best.auc.list, 
                      steatosis.bestmodels$best.auc.list)

saveRDS("best.auc.list")

# Combine results into a final dataframe
best.auc.df <- rbind(F1.bestmodels$best.auc.df, 
                     F.bestmodels$best.auc.df, 
                     Inflammation.bestmodels$best.auc.df, 
                     steatosis.bestmodels$best.auc.df)

saveRDS("best.auc.df")

# Summarize best models features 
# ---------------------------------------------------------------------------

F1.topFeatures<-reshape2::melt(F1.bestmodels$topFeatures.list)
colnames(F1.topFeatures)<-c("Feature","n.Features", "omic")
F1.topFeatures$panel<-"F>=2"

F.topFeatures<-reshape2::melt(F.bestmodels$topFeatures.list)
colnames(F.topFeatures)<-c("Feature","n.Features", "omic")
F.topFeatures$panel<-"F>=3"

Inflammation.topFeatures<-reshape2::melt(Inflammation.bestmodels$topFeatures.list)
colnames(Inflammation.topFeatures)<-c("Feature","n.Features", "omic")
Inflammation.topFeatures$panel<-"I>=2"

Steatosis.topFeatures<-reshape2::melt(steatosis.bestmodels$topFeatures.list)
colnames(Steatosis.topFeatures)<-c("Feature","n.Features", "omic")
Steatosis.topFeatures$panel<-"S>=2"

topFeatures<-rbind(F1.topFeatures,F.topFeatures,Inflammation.topFeatures,Steatosis.topFeatures)

# merge with best model table
topFeatures$omic.model<-paste(topFeatures$omic, topFeatures$n.Features,topFeatures$panel, sep = ".")
best.auc.df$omic.model<-paste(best.auc.df$omic,"top" ,best.auc.df$no.features, best.auc.df$panel,sep=".")

topFeatures.bestmodel<- topFeatures %>% dplyr::filter(omic.model %in% best.auc.df$omic.model)
topFeatures.bestmodel<-base::merge(topFeatures,best.auc.df, by="omic.model")

topFeatures.bestmodel <- topFeatures.bestmodel %>% dplyr:::select.data.frame(omic.model,omic.x,n.Features,Feature,AUC.cv, sd,panel.x)
topFeatures.bestmodel$omic.name<-str_split_fixed(topFeatures.bestmodel$omic.x,"\\.",2)[,1]

#saveRDS(topFeatures.bestmodel,"/projects/arumugam/people/xhf865/Camila/GALAXY/GALA/XGboost_2025/Resuts/Best.Models/topFeatures.bestmodel.plots.rds")

topFeatures.bestmodel$Feature<-str_split_fixed(topFeatures.bestmodel$Feature, "\\.",2)[,2]

topFeatures.bestmodel.summary<- topFeatures.bestmodel %>% 
  group_by(omic.model) %>% dplyr::summarize(omic.name,     
                                            panel=panel.x,
                                            n.features=n.Features,
                                            Features= paste0(Feature, collapse = ', '), 
                                            AUC.cv=AUC.cv,
                                            sd=sd)

topFeatures.bestmodel.summary<-unique(topFeatures.bestmodel.summary)

# saveRDS(topFeatures.bestmodel.summary)





