# Author: Camila Alvarez-Silva

library(dplyr)
library(ggplot2)
library(ggrepel)
library(stringr)

# Function: lin_reg_associate
# This function performs an association analysis between omics data and clinical outcomes,
# correcting for specified confounders. It applies linear regression models, corrects for multiple testing,
# and provides summary statistics along with visualization plots.

lin_reg_associate <- function(df.matrix, metadata, feature, confounders) {
  
  # Check for missing inputs
  if (missing(df.matrix) || missing(metadata) || missing(feature) || missing(confounders)) {
    stop("Missing required input")
  }
  
  # Extract unique omic types from column names
  omics <- unique(str_split_fixed(colnames(df.matrix), "\\.", 2)[,1])
  
  # Define color palette for plots
  color_levels <- c("#175e63", "#52796F", "#E7B800", "#FC4E07", "#ed070b", "#2A9D8F", "#74a94c",
                    "#F4A261", "#F3DFC1", "#46ACAF", "#143F52", "#247294", "#A59132", "#CAD2C5",
                    "#84A98C", "#5c7b99", "#396a68", "#E76F51")
  names(color_levels) <- omics
  
  df.list <- list()
  
  # Iterate over each omic type
  for (name in omics) {
    
    # Extract relevant columns from the omics dataset
    name.df <- df.matrix[, startsWith(colnames(df.matrix), name)]
    
    # Remove rows where all values are NA or zero
    if (any(rowSums(is.na(name.df)) == ncol(name.df))) {
      name.df <- name.df[rowSums(is.na(name.df)) != ncol(name.df), ]
    }
    
    name.df$Sample_ID <- rownames(name.df)
    
    # Prepare clinical feature dataframe
    feature.df <- metadata %>% dplyr::select(SampleID, all_of(confounders), all_of(feature)) %>% 
      filter(!is.na(.data[[feature]]))
    colnames(feature.df)[colnames(feature.df) == "SampleID"] <- "Sample_ID"
    colnames(feature.df)[colnames(feature.df) == feature] <- "feature"
    
    # Merge omics and feature dataset
    name.df.feature <- merge(feature.df, name.df, by = "Sample_ID")
    rownames(name.df.feature) <- name.df.feature$Sample_ID
    name.df.feature$Sample_ID <- NULL
    
    # Remove columns where all values are zero or NA
    if (any(colSums(name.df.feature, na.rm = TRUE) == 0)) {
      name.df.feature <- name.df.feature[, colSums(name.df.feature, na.rm = TRUE) != 0]
    }
    
    # Standardize data
    name.df.feature.scale <- as.data.frame(scale(name.df.feature))
    rownames(name.df.feature.scale) <- rownames(name.df.feature)
    
    # Apply linear models
    lm.out <- apply(name.df.feature.scale[, -(1:(length(confounders)+1))], 2, function(x) {
      ## Create a data frame including confounders and clinical feature of interest
      df <- data.frame(Feature = x, name.df.feature.scale[, 1:(length(confounders)+1)]) 
      summary(lm(Feature ~ ., data = df))
    })
    
    # Extract coefficient statistics for the feature of interest
    lm.out <- lapply(lm.out, function(x) x$coefficients["feature", ])
    lm.df <- as.data.frame(do.call("rbind", lm.out))
    colnames(lm.df)[4] <- "p.val"
    lm.df$omic <- name
    lm.df$feature<-str_split_fixed( row.names(lm.df) , "\\.", 2)[,2]
    
    df.list[[name]] <- lm.df
  }
  
  # Combine results and adjust for multiple testing
  df.feature.omics <- do.call("rbind", df.list)
  row.names(df.feature.omics)<-str_split_fixed(row.names(df.feature.omics), "\\.",2)[,2]
  df.feature.omics$p.adjust <- p.adjust(df.feature.omics$p.val, method = "fdr")
  
  # Filter significant results
  df.feature.omics.filt <- df.feature.omics %>% filter(p.adjust < 0.05)
  df.feature.omics.filt$outcome <- feature
  
  # Generate p-value plot
  mean.pval <- df.feature.omics.filt %>% group_by(omic) %>% 
    summarise(mean.pval = mean(p.adjust, na.rm = TRUE)) %>% 
    arrange(mean.pval)
  mean.pval$omic <- factor(mean.pval$omic, levels = mean.pval$omic)
  df.feature.omics.filt$omic <- factor(df.feature.omics.filt$omic, levels = rev(mean.pval$omic))
  
  plot.pval <- ggplot(df.feature.omics.filt, aes(x = -log10(p.adjust), y = omic, fill = omic)) +
    geom_point(shape = 21, alpha = 0.65) +
    theme_bw() +
    theme(legend.position = "none", text = element_text(size = 11)) +
    geom_vline(xintercept = -log10(0.05), col = "black") +
    xlab("-log10(p.adjust)") + ylab(NULL) +
    facet_grid(~outcome) +
    scale_fill_manual(values = color_levels)
  
  # Generate volcano plot
  df.feature.omics$label <- ifelse(df.feature.omics$p.adjust < 0.05, df.feature.omics$feature, NA)
  
  
  # Select top 2 features per omic (top positive and negative associoations)
  label.list<-list()
  omics<-unique(df.feature.omics.filt$omic)
  df.feature.omics$omic.feature<-paste(df.feature.omics$omic,df.feature.omics$feature, sep=".")
  
  for(i in omics){
    
    df.filt<-filter(df.feature.omics,omic ==i)
    
    #Take the top  features by effect size 
    df.filt.neg<-df.filt %>% filter(p.adjust < 0.05 & Estimate < 0)
    df.filt.neg<-df.filt.neg[order(abs(df.filt.neg$Estimate),decreasing = T),]
    
    
    df.filt.pos<-df.filt %>% filter(p.adjust < 0.05 & Estimate > 0)
    df.filt.pos<-df.filt.pos[order(df.filt.pos$Estimate,decreasing = T),]
    
    label.list[[i]]<-c(df.filt.pos$omic.feature[1:2],df.filt.neg$omic.feature[1:2])
  }
  
  label.df<-reshape2::melt(label.list) %>% filter(!is.na(value))
  label<-label.df$value
  
  
  df.feature.omics$label<-ifelse(df.feature.omics$omic.feature %in% label, df.feature.omics$feature,NA)
  
  # Edit label name
  df.feature.omics$label<-str_remove_all(df.feature.omics$label,"_Relative Amount")
  df.feature.omics$label<-str_remove_all(df.feature.omics$label," incertae sedis")
  df.feature.omics$label<-str_split_fixed(df.feature.omics$label, "\\[meta_",2)[,1]
  df.feature.omics$label<-str_split_fixed(df.feature.omics$label, "\\[ref_",2)[,1]
  df.feature.omics$label<-str_split_fixed(df.feature.omics$label, "\\: ",2)[,1]
  
  
  
  plot.volcanoPlot <- ggplot(df.feature.omics, aes(x = Estimate, y = -log10(p.adjust), fill = omic)) +
    geom_point(shape=21, alpha=0.65) + 
    geom_hline(yintercept=-log10(0.05), col="black",linetype="dotted")+
    geom_label_repel(aes(label = label, fill = omic), 
                     color = "white", 
                     fontface = "bold", 
                     box.padding = 0.3, 
                     point.padding = 0.2, 
                     segment.color = "black", 
                     max.overlaps = Inf,
                     size = 2.5) +
    scale_fill_manual(values = color_levels) +
    theme_bw() +
    xlab("Standardized beta coefficient")+
    theme(legend.position = "right",
          panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(),
          panel.background = element_blank(),
          axis.text=element_text(size=9),
          axis.ticks = element_blank(),
          text = element_text(size = 12))
  
  
  
  # Return results
  return(list(df.feature.omics = df.feature.omics, 
              plot.pval = plot.pval, 
              plot.volcanoPlot = plot.volcanoPlot))
}
