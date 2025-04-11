## R script to construct prognostic model for clinical outcomes in the GALA-ALD
## Author: Suguru Nishijima

## Load library
library(tidyverse)
library(mlr)
library(survival)
library(glmnet)
library(rsample)
library(data.table)

## Function for 5-fold cross-validation and predictive scores
my_cv <- function(d, surv, data, type){
  set.seed(1)
  res <- list()

  d2 <- d
  colnames(d2) <- colnames(d2) %>% make.names(unique = T)
  
  ## Add surv data to d
  temp.day <- surv %>% as.character() %>% str_replace("\\+", "") %>% as.character() %>% as.numeric()
  temp.cond <- ifelse(grepl("\\+", as.character(surv)), F, T)
  d2 <- data.frame(d2, day = temp.day, cond = temp.cond)
  
  ## 5 times repeated 5-fold cross-validation
  sp.df <- c()
  pred.df <- c()
  feat.list <- list()
  cindex.feat.df <- c()
  all_df <- c()
  
  j <- 1
  for(j in 1:5){
    set.seed(j)
    train_folds <- vfold_cv(d2, v = 5)
    
    ## 5-fold cross-validation
    feat.list.temp <- list()  
    i <- 1
    for(i in 1:5){
      print(paste0("rep: ", j, ", fold: ", i))
      
      ## Prepare training data (80% of original dataset)
      x1 <- analysis(train_folds$splits[[i]])
      
      ## Prepare test data (20% of original dataset)
      x2 <- assessment(train_folds$splits[[i]])
      colnames(x1) %>% sort() %>% tail()
      
      ## Define model
      task <- makeSurvTask(data = x1, target = c("day", "cond"))
      lrn.ranger <- makeLearner("surv.ranger", importance = "permutation")

      ## Construct model
      model <- mlr::train(lrn.ranger, task)
      model$features %>% sort() %>% tail()
      
      ## Feature selection
      ## Calculate feature importance
      temp <- model %>% mlr::getFeatureImportance()
      model.imp <- temp$res %>% data.frame()
      colnames(model.imp) <- c("rowname", "imp")
      model.imp <- model.imp %>% arrange(-imp)
        
      if(ncol(d) < 15){
        feat_num_max <- ncol(d)
      }else{
        feat_num_max <- 15
      }

      cindex.vec <- c()
      pred.df.fs <- c()      
      feat_num <- 2
      for(feat_num in 2:feat_num_max){
        paste0("fs: ", feat_num) %>% print()
        
        ## Extract only top features and construct model
        x1.red <- x1[, model.imp$rowname[1:feat_num]] %>% mutate(day = x1$day) %>% mutate(cond = x1$cond)
        task = makeSurvTask(data = x1.red, target = c("day", "cond"))
        lrn.ranger.red = makeLearner("surv.ranger")
        model.red = mlr::train(lrn.ranger.red, task)
        
        ## Apply the model to unused data
        x2.surv <- Surv(x2$day, x2$cond)
        pred <- predict(model.red, newdata = x2)
        pred <- pred$data %>% data.frame(surv = x2.surv, rep = j, fold = i, fs = feat_num) %>% rownames_to_column()
        pred <- pred %>% select(-contains("truth"))
        pred.df.fs <- rbind(pred.df.fs, pred) %>% arrange(rowname)
        
        ## Evaluate the performance using Cindex
        cindex.vec[feat_num] <- Cindex(pred$response, x2.surv)
      }

      cindex.feat.df.temp <- data.frame(cindex = cindex.vec, feat_num = 1:feat_num_max, i = i, j = j)
      cindex.feat.df <- rbind(cindex.feat.df, cindex.feat.df.temp)
      feat.list.temp[[i]] <- model.imp$rowname
      ##
      
      ## Apply the  model to test data
      pred <- predict(model, newdata = x2)
      pred <- pred$data %>% data.frame(surv = x2.surv, rep = j, fold = i, fs = "all") %>% rownames_to_column()
      pred.df <- pred %>% select(-contains("truth"))
      
      ## Save results for all and selected features
      all_df <- rbind(all_df, pred.df, pred.df.fs)
    }
    feat.list[[j]] <- feat.list.temp
  }
  
  ## Take average across cross-validations
  pred.mean <- tapply(all_df$response, all_df$rowname, mean)
  pred.median <- tapply(all_df$response, all_df$rowname, median)
  
  ## Save results
  lrn.ranger = makeLearner("surv.ranger")
  task = makeSurvTask(data = d2, target = c("day", "cond"))
  fullmodel = mlr::train(lrn.ranger, task)
  
  res[[1]] <- pred.mean
  res[[2]] <- glmnet::Cindex(pred.mean, surv)
  res[[3]] <- fullmodel
  res[[4]] <- surv
  res[[5]] <- d
  res[[6]] <- d2  
  res[[7]] <- c
  res[[8]] <- cindex.feat.df
  res[[9]] <- feat.list
  res[[10]] <- all_df
  
  output <- paste0("out/prognostic_model/", data, ".", type, ".rds")
  saveRDS(res, file = output)
  return(res)
}



## Read outcome data
md <- readRDS("data/GALA-ALD.clinical_outcomes.rds")

## Read omic data to analyze
files <- list.files("data/omic_data/") %>% paste0("data/omic_data/", .)
files

## Add clinical data
files <- c(files, "data/metadata/GALA_ALD.clinic.tsv", "data/metadata/GALA_ALD.confounder.tsv")
files


## For loop for all files
res.df <- c()
for(i in 1:length(files)){
  
  file <- files[i]
  print(i)  
  print(file)

  data <- file %>% str_replace("..*\\/", "") %>% str_replace("\\..*", "")
  
  ## Read file
  d <- fread(file) %>% as.matrix(rownames = 1) %>% data.frame(check.names = F)
  
  if(grepl("clinic", file)){
    data <- "clinic"
  }
  if(grepl("confounder", file)){
    data <- "confounder"
  }

  ## Exclude feature with only NA
  keep <- apply(is.na(d), 2, sum) != nrow(d)
  d <- d[, keep]

  ## Extract samples with outcome data (GALA-ALD)
  keep1 <- which(rownames(d) %in% md$Sample_ID)
  keep2 <- md$Sample_ID %in% rownames(d)
  d2 <- d[keep1, ]
  md2 <- md[keep2, ]
    
  ## Replace NA with half of lowest value
  if(sum(is.na(d2)) > 0){
    for(ii in 1:ncol(d2)){
      na.pos <- is.na(d2[, ii])
      d2[na.pos, ii] <- (min(d2[, ii], na.rm = T) / 2)
    }
  }
  
  ## Exclude Shannon diversity
  shannon.sum <- grepl("Shannon", colnames(d2)) %>% sum()
  if(shannon.sum > 0){
    d2 <- d2 %>% select(-Shannon)    
  }
  
  if(nrow(d2) == 0){next}
  
  ## Construct prognostic model
  surv <- md2$decomp_surv
  d <- d2
  
  res.liver <- my_cv(d2, md2$decomp_surv, data, "decompensation")
  res.death <- my_cv(d2, md2$death_surv, data, "death")
  res.inf <- my_cv(d2, md2$inf_surv, data, "infection")
}
