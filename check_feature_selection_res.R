## R script to check results of feature selection
## Author: Suguru Nishijima

## Load library
library("tidyverse")

## Read output files from construct_prognostic_model.R
files <- list.files("out/prognostic_model/", pattern = "rds$") %>% paste0("out/prognostic_model/", .)
files <- files[!grepl("evaluation", files)]

df <- c()
df.feat <- c()
df.error <- c()
for(file in files){
  print(file)
  omics <- file %>% str_replace(".*\\/", "") %>% str_replace("\\..*", "")
  temp <- file %>% str_split("\\.", simplify = T)
  outcome <- temp[2]
  res <- read_rds(file)
  
  ## C-index in each iteration of feature selection
  temp.res <- res$cindex_per_featnum %>% group_by(feat_num) %>% summarise(m = mean(cindex, na.rm = T)) %>% arrange(-m) %>% data.frame()
  m.max <- temp.res$m %>% max(na.rm = T)
  
  ## Identify highest number of features keeping 99% of cindex 
  m.base <- temp.res$m[1]
  feat_num_max <- temp.res$feat_num[1]
  best.num <- 2
  for(i in feat_num_max:2){
    print(i)
    if(m.base * 0.99 > temp.res$m[temp.res$feat_num == i]){
      best.num <- i + 1
      break
    }
  }
  
  ## Check selected features
  temp <- res$cindex_per_featnum %>% filter(feat_num == best.num) %>% arrange(-cindex)
  best.i <- temp$i[1]
  best.j <- temp$j[1]
  
  feat.list <- res$feature_list[[best.i]][[best.j]]
  feat.list[1:best.num] %>% print()

  cindex <- temp.res$m[1] ## mean cindex of the best num
  res.max <- res$cindex_per_featnum %>% filter(feat_num == best.num) %>% arrange(-cindex) %>% head(1)
  
  feat.list <- res$feature_list[[res.max$i]][[res.max$j]]
  feat.list <- feat.list[1:res.max$feat_num]

  ## Save c-index with the best features  
  temp <- data.frame(omics = omics, outcome = outcome, cindex = cindex, best_num = best.num)
  df <- rbind(df, temp)
  
  ## Save best features
  temp <- data.frame(omics = omics, feat = feat.list, best_num = best.num, outcome = outcome, i = res.max$i, j = res.max$j)
  df.feat <- rbind(df.feat, temp)
  
  ## Save c-index of each iteration with the best features
  c.error <- res$cindex_per_featnum %>% filter(feat_num == best.num) %>% group_by(feat_num, i) %>% summarise(m = mean(cindex))
  temp <- data.frame(omics = omics, cindex = c.error$m, best_num = best.num, outcome = outcome)  
  df.error <- rbind(df.error, temp)
}

## Save results
saveRDS(df.feat, file = "out/rds/selected_feature.rds")
saveRDS(df, file = "out/rds/df.selected_feature.rds")

