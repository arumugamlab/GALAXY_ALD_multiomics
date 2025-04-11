## R script for cox regression analysis of the GALA-ALD samples
## Author: Suguru Nishijima

## load library
library(tidyverse)
library(survival)
library(data.table)

## read metadata
md <- readRDS("data/GALA-ALD.clinical_outcomes.rds")
md.conf <- readRDS("data/metadata/GALA_ALD.confounder.rds")
md.clinic <- readRDS("data/metadata/GALA_ALD.clinic.rds")
md.marker <- readRDS("data/metadata/GALA_ALD.marker.rds")

## use only shared samples
keep <- rownames(md.conf) %in% md$Sample_ID
md.conf <- md.conf[keep, ]
md.clinic <- md.clinic[keep, ]
md.marker <- md.marker[keep, ]

## read omic data to analyze
files <- list.files("data/omic_data/") %>% paste0("data/omic_data/", .)
files

res.df <- c()
my_cox <- function(num, md, data){
  file <- files[num]
  print(num)  
  print(file)    
  d <- fread(file) %>% as.matrix(rownames = 1) %>% data.frame(check.names = F)
  
  ## extract only shared ALD data between omic and surv datasets
  keep <- which(rownames(d) %in% md$Sample_ID)
  d <- d[keep, ]
  
  keep <- which(md$Sample_ID %in% rownames(d))
  if(length(keep) == 0){return()}  
  md2 <- md[keep, ]
  md.conf2 <- md.conf[keep, ]
  
  ## select outcome data to use
  if(grepl("decompensation", data)){surv <- md2$decomp_surv}
  if(grepl("death", data)){surv <- md2$death_surv}
  if(grepl("inf", data)){surv <- md2$inf_surv}  
  
  ## for loop for cox regression
  hz <- c()
  p.cox <- c()  
  for(i in 1:ncol(d)){
    
    if(sum(is.na(d[, i])) == nrow(d)){
      hz[i] <- NA
      p.cox[i] <- NA
      next
    }
    
    ## Cox regression with adjustment for confounders
    d.temp <- data.frame(surv = surv, feat = d[, i], md.conf2, check.names = F) %>% data.frame()
    cox.res <- coxph(surv ~ ., data = d.temp, iter.max = 50) %>% summary()

    hz[i] <- cox.res$conf.int["feat", "exp(coef)"]
    p.cox[i] <- cox.res$coefficients["feat", "Pr(>|z|)"]  
  }
  
  ## save results
  res <- data.frame(name = colnames(d), hz = hz, p = p.cox, outcome = data) %>% arrange(p)
  file2 <- file %>% str_replace(".*\\/", "")
  
  output <- paste0("out/cox_regression/", file2, ".", data, ".cox.tsv")
  write.table(res, file = output, sep = "\t", col.names = NA)
  
  return(res)
}

## for all the omics datasets
for(i in 1:length(files)){
  print(files[i])
  res.liver <- my_cox(i, md, "decompensation")
  res.death <- my_cox(i, md, "death")
  res.surv <- my_cox(i, md, "infection")
}

