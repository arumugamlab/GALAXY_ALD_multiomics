## R script to evaluate prognostic models constructed
## Author: Suguru Nishijima

## load library
library(tidyverse)
library(survival)

## Function to compute various performance metrics
calc_various_metrics <- function(d){
  time <- 365 * 5
  nri.cut <- 0.5
  
  ## c-index
  c <- glmnet::Cindex(d$marker, d$outcome)
  
  ## time-dependent AUC
  res.roc <- timeROC::timeROC(T = d$days, delta = d$status, marker = d$marker, cause = 1, times = c(365, 365*3, 365*5))
  roc.1yr <- res.roc$AUC[1]
  roc.3yr <- res.roc$AUC[2]
  roc.5yr <- res.roc$AUC[3]
  
  ## Net Reclassification Index (NRI)
  d2 <- d[!is.na(d$te), ]
  
  ## Fit models using marker and TE separately
  mod1 <- coxph(d2$outcome ~ d2$marker, x = T)
  fit.marker <- survfit(mod1, d2)
  d2$marker.prob <- summary(fit.marker, times = time)$surv %>% as.numeric()
  
  mod2 <- coxph(d2$outcome ~ d2$te, x = T)
  fit.te <- survfit(mod2, d2)
  d2$te.prob <- summary(fit.te, times = time)$surv %>% as.numeric()
  
  ## Run NRI calculation
  sink(tempfile())
  suppressMessages({
    temp <- nricens::nricens(time = d2$days, event = d2$status, mdl.std = mod2, mdl.new = mod1, updown = "diff", cut = 0, t0 = 365, niter = 0, msg = F)    
    nri.1yr <- temp$nri$Estimate[1] 
    
    temp <- nricens::nricens(time = d2$days, event = d2$status, mdl.std = mod2, mdl.new = mod1, updown = "diff", cut = 0, t0 = 365*3, niter = 0, msg = F)    
    nri.3yr <- temp$nri$Estimate[1] 
    
    temp <- nricens::nricens(time = d2$days, event = d2$status, mdl.std = mod2, mdl.new = mod1, updown = "diff", cut = 0, t0 = 365*5, niter = 0, msg = F)    
    nri.5yr <- temp$nri$Estimate[1] 
  })
  sink()
  
  ## summarise results
  res <- data.frame(Cindex = c, 
                    
                    AUC.1yr = roc.1yr,
                    AUC.3yr = roc.3yr,
                    AUC.5yr = roc.5yr, 
                    
                    NRI.1yr = nri.1yr,
                    NRI.3yr = nri.3yr,
                    NRI.5yr = nri.5yr
  )
  res <- res %>% mutate(n_total = nrow(d2))
  
  return(res)
}


## Read output files from construct_prognostic_model.R
files <- list.files("out/prognostic_model/", pattern = "rds$") %>% paste0("out/prognostic_model/", .)
files

## Read clinical marker data
md <- read_rds("data/metadata/GALA_ALD.marker.rds")
md$sampleID <- rownames(md)

# -----------------------------------------
# Section 1: Evaluate full models (using all features)
# -----------------------------------------

res_all <- c()
for(file in files){
  print(file)
  
  ## Extract omic type and outcome type from file name
  omics <- file %>% str_replace(".*\\/", "") %>% str_replace("\\..*", "")
  temp <- file %>% str_split("\\.", simplify = T)
  outcome <- temp[2]
  
  ## Load model results
  model_res <- read_rds(file)
  
  ## Extract predictions for full models (i.e., all features)
  df <- model_res$prediction_df
  df <- df %>% filter(fs == "all")

  res <- c()
  ## Loop through each of the 5 cross-validation repetitions
  for(r in 1:5){
    paste0(omics, ", ", outcome, ", rep: ", r) %>% print()
    
    ## for each rep
    d <- df %>% filter(rep == r) %>% arrange(rowname)
    
    surv <- Surv(d$surv[, 1], d$surv[, 2])
    
    md.sub <- data.frame(sampleID = d$rowname, marker = d$response, outcome = surv)
    md.sub <- left_join(md.sub, md, by = "sampleID")
    
    md.sub$days <- md.sub$outcome %>% str_remove("\\+") %>% as.numeric()
    md.sub$status <- ifelse(grepl("\\+", md.sub$outcome), F, T)
    d <- md.sub
    
    ## Calculate performance metrics
    temp <- calc_various_metrics(md.sub)
    
    ## Store results with annotation
    temp <- temp %>% data.frame(., omics = omics, outcome = outcome, rep = r, fold = NA)
    res <- rbind(res, temp)
  }
  res_all <- rbind(res_all, res)
}
df <- res_all

## Save results
rownames(df) <- c()
df$outcome <- df$outcome %>% 
  str_replace("decompensation", "Decompensation") %>% 
  str_replace("death", "All-cause mortality") %>% 
  str_replace("infection", "Infection")

saveRDS(df, file = "out/evaluation_metrics.full_model_with_each_repetition.rds")
write_tsv(df, file = "out/evaluation_metrics.full_model_with_each_repetition.tsv")


# -----------------------------------------
# Section 2: Evaluate reduced models (with feature selection)
# -----------------------------------------


## summarize results for best feature-selected model with each repetition
## read selected features
d.best_feat <- read_rds("out/rds/selected_feature.rds") %>% select(omics, best_num, outcome, i, j) %>% unique()
d.best_feat

res_all <- c()
for(file in files){
  print(file)
  
  ## Extract omic type and outcome type from file name
  file.omics <- file %>% str_replace(".*\\/", "") %>% str_replace("\\..*", "")
  temp <- file %>% str_split("\\.", simplify = T)
  file.outcome <- temp[2]
  
  ## Load model results
  model_res <- read_rds(file)

  ## Get the best number of features selected for this omic-outcome combination
  d.best <- d.best_feat %>% filter(omics == file.omics & outcome == file.outcome)
  
  df <- model_res$prediction_df
  df <- df %>% filter(fs == d.best$best_num)
  
  f <- 1
  res <- c()
  for(r in 1:5){
    paste0(file.omics, ", ", file.outcome, ", rep: ", r, ", feature_num: ", d.best$best_num) %>% print()
    d <- df %>% filter(rep == r) %>% arrange(rowname)
    surv <- Surv(d$surv[, 1], d$surv[, 2])
    
    md.sub <- data.frame(sampleID = d$rowname, marker = d$response, outcome = surv)
    md.sub <- left_join(md.sub, md, by = "sampleID")
    
    md.sub$days <- md.sub$outcome %>% str_remove("\\+") %>% as.numeric()
    md.sub$status <- ifelse(grepl("\\+", md.sub$outcome), F, T)
    d <- md.sub
    
    ## Calculate performance metrics
    temp <- calc_various_metrics(md.sub)
    
    ## Store results with annotation
    temp <- temp %>% data.frame(., omics = file.omics, outcome = file.outcome, rep = r, fold = NA, best_num = d.best$best_num)
    res <- rbind(res, temp)
  }
  res_all <- rbind(res_all, res)
}
df <- res_all

## Save results
rownames(df) <- c()
df$outcome <- df$outcome %>% 
  str_replace("decompensation", "Decompensation") %>% 
  str_replace("death", "All-cause mortality") %>% 
  str_replace("infection", "Infection")

saveRDS(df, file = "out/evaluation_metrics.reduced_model_with_each_repetition.rds")
write_tsv(df, file = "out/evaluation_metrics.reduced_model_with_each_repetition.tsv")
