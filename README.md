# 1. Phenotype-omics associations using `lin_reg_associate()`

### Overview

This R function `lin_reg_associate()` perform association analysis between omics data and clinical outputs, correcting by given confounders. The function applies linear regression models and corrects for multiple testing, providing summary statistics and visualization plots.
Function used to generate Figure 1.

### Dependencies

Ensure these packages are installed:

```r
install.packages(c("dplyr", "ggplot2", "ggrepel", "stringr"))
```

### Inputs

| Parameter                          | Type         | Description |
|---------------------------------|-----------------|-------------|
| `df.matrix`                     | DataFrame     | Data matrix with omics features (columns) and samples (rows). **Note:** each feature must be prefixed with its corresponding omic type (e.g.,"Proteomics.P10643") |
| `metadata`                      | DataFrame     | Contains clinical information related to the samples. **Note:** make sure column with sample ID is named as "SampleID"|
| `Feature`                       | String        |The column name in the metadata representing the clinical feature of interest. Note: The feature must be a numeric variable (e.g., "steatosis_numeric").|
| `confounders`                   | Vector (String)   |A vector containing the column names in the metadata for the variables to be used as confounders.|

### Outputs

| Outputs                         | Type            | Description |
|---------------------------------|-----------------|-------------|
| `df.feature.omics`              | DataFrame    | Data frame with the linear model association results between omics features and the  clinical feature of interest. |
| `plot.pval`                     | ggplot object   | A scatter plot displaying the p-adjusted values of the associations stratified by omic.|
| `plot.volcanoPlot`              | ggplot object   |A volcano plot visualizing the effect size and significance of associations, colored by omic|

### Example Usage
```r

Result.associations <- lin_reg_associate(
    df.matrix = df.matrix,
    metadata = metadata,
    feature="steatosis_numeric",
    confounders=c("Age","BMI","gender")
    )

# Display results
print(Result.associations$df.feature.omics)
print(Result.associations$plot.pval)
print(Result.associations$plot.volcanoPlot)

```
# 2. Biomarker discovery

We use the [XGBoost](https://doi.org/10.1145/2939672.2939785) machine learning algorithm to derive compact biomarker sets for a given phenotype. We typically set aside a test set, and tune hyperparameters using `n` times `n`-fold cross-validation on the training set. Once the final hyperparameter set is selected, the optimal model is then evaluated on the test set -- once and only once. Here are the `R` functions that help you with that.

## Training XGBoost model with `xgboost_train()`

### Overview

This R function implements an XGBoost-based classification pipeline, including feature reduction and performance evaluation via cross-validation. It supports multiple omics datasets and outputs a range of performance metrics and model details.

### Dependencies

Ensure these packages are installed:

```r
install.packages(c("dplyr", "xgboost", "caret", "pROC", "Metrics", "reshape2", "ggplot2", "xpectr"))
```

### Inputs

| Parameter        | Type              | Description |
|-----------------|------------------|-------------|
| `omics`         | Character Vector  | Types of omics data (e.g., "Metabolomics", "Metagenomics", "Proteomics"). |
| `feature`       | DataFrame         | Target labels for classification. |
| `seed`          | Integer           | Random seed for reproducibility. |
| `top`           | Integer           | Number of top features to retain in feature reduction. |
| `df.matrix`     | DataFrame         | Data matrix with omics features (columns) and samples (rows). **Note:** each feature must be prefixed with its corresponding omic type (e.g.,"Proteomics.P10643")|
| `panel`         | String            | Column name on the metadata of the feature that is going to be analyzed (e.g., "PreACLF"). |
| `nfolds`        | Integer           | Number of folds for cross-validation. |
| `ntimes`        | Integer           | The amount of repeated fold computations to perform. |
| `Microbiome.features` | Character Vector  | List of microbiome-related features. |
| `Host.features` | Character Vector  | List of host-related features. |
| `metadata`      | DataFrame         | Metadata containing sample IDs and feature of interest. |
| `standards.variables` | Character Vector |Standards variables names to be evaluated individually. (e.g., "TE", "Meld.score", "AST")|

### Outputs

The function returns a list of outputs covering model performance, predictions, and feature importance:

| Output                         | Type            | Description |
|--------------------------------|----------------|-------------|
| `loss.list`                    | List (Plots)   | Loss curves for each omic type, showing training vs validation loss over one fold. |
| `cvAUC.list`                   | List           | Cross-validated AUC scores for each omic type. |
| `cvAUC.list.top`               | List           | Cross-validated AUC scores using only the top X most important features. |
| `feature.importance.list`       | List (DataFrames) | Feature importance for each fold and omic type from the full dataset model. |
| `feature.importance.list.top`   | List (DataFrames) | Feature importance using only the reduced top features. |
| `cvAUC.value`                  | List (Numeric) | Final AUC values after cross-validation for each omic type. |
| `cvAUC.se`                     | List (Numeric) | Standard error of AUC across cross-validation folds. |
| `model_list.all`                | List (Models)  | Trained XGBoost models for each fold and omic type. |
| `model_list.all.top`            | List (Models)  | XGBoost models trained using only the top features. |
| `predictions.list`              | List (Vectors) | Predicted probabilities for each fold’s validation set. |
| `predictions.list.top`          | List (Vectors) | Predicted probabilities from models trained on top features only. |
| `labelsValidation.list`         | List (Vectors) | True validation set labels for each fold. |
| `folds`                         | List (Vectors) | Indices of training/validation samples for each cross-validation fold. |
| `out`                           | List           | AUC and confidence intervals for the complete model. |
| `out.top`                       | List           | AUC and confidence intervals for the top feature model. |
| `score_list.all`                | List (DataFrames) | AUC scores, rounds, and evaluation logs per fold. |
| `confusionMatrix.list.all`      | List (Matrices) | Confusion matrices for each fold’s predictions. |
| `confusionMatrix.lit.top.all`   | List (Matrices) | Confusion matrices from top feature models. |
| `train.labels.ID.out`           | List (Vectors) | Training set sample IDs per fold. |
| `val.labels.ID.out`             | List (Vectors) | Validation set sample IDs per fold. |

### Example Usage

```r
result <- xgboost_train(
  omics =  c("Proteomics","Metabolomics","Microbiome", "Host"),
  feature = feature_df,
  seed = 123,
  top = 15,
  df.matrix = df_matrix,
  panel = "PreACLF",
  nfolds = 5,
  ntime = 5,
  Microbiome.features = Microbiome.features,
  Host.features = Host.features,
  metadata = metadata_df
)

# Access results:
result$cvAUC.list          # Full model AUCs
result$cvAUC.list.top      # Top feature model AUCs
result$feature.importance.list # Feature importance for full models
result$feature.importance.list.top # Feature importance for top models
```


## Evaluating a trained XGBoost model with `xgboost_eval()`

### Overview

This R function  `xgboost_eval()` generates performance plots for a trained XGBoost model and outputs relevant performance metrics and visualizations.


### Inputs

| Parameter       | Type                  | Description |
|-----------------|-----------------------|-------------|
| `xgboost.object` | Output of `xgboost_train()` | A trained XGBoost model with AUC values and feature importance information. |
| `top`           | Integer               | Number of top features to include in visualizations. |

### Outputs

| Output                          | Type            | Description |
|---------------------------------|-----------------|-------------|
| `plot.complete`                 | List (Plots)    | Complete set of plots combining AUC curves and feature upset plots. |
| `AUC.top.feature.curve`         | Plot            | Line plot showing AUC scores for top features across models. |
| `plot.AUC.all.feature.reduced`  | List (Plots)    | AUC bar plots for selected top features. |
| `AUC.barplot.all.omics`         | Plot            | Bar plot summarizing AUC values for all omics data. |
| `confusionMatrix.best.model`    | List (Plots)    | Confusion matrix for the best-performing model across feature reduced models. |
| `bestModel`                     | DataFrame       | Summary table of the best-performing models per omic. |
| `topFeatures.list`              | List (Vectors)  | List of top features from the best models. |
| `AUC.barplot.omics.all.bestmodels` | Plot         | AUC bar plot for best models across all omics. |
| `bestModel.top`                 | DataFrame       | Table summarizing the minimum number of features in the best-performing models. |
| `AUC.all.top`                   | DataFrame       | Melted data frame of AUC values across top models and features. |
| `cvAUC.list.melt`               | DataFrame       | Refined data frame of AUC values with features reordered for visualization. |

### Summary Code Breakdown

The function performs the following steps:

1. **Generate AUC Bar Plots:**
   - Calculates cross-validation AUC values and their standard deviations.
   - Creates bar plots for each omic, displaying AUC scores and error bars.

2. **Plot Top Features:**
   - Extracts performance metrics for models with varying top features.
   - Generates AUC plots with error bars and labels showing top feature performance.

3. **Best Model Identification:**
   - Identifies the best model per omic based on maximum AUC and minimun number of features.
   - Extracts feature importance tables for each best-performing model.

4. **Confusion Matrix Visualization:**
   - Creates confusion matrix plots with color-coded correct/incorrect predictions.

5. **Upset Plot Generation:**
   - Uses the `ComplexUpset` package to generate upset plots for feature overlaps across feature reduced models.

### Dependencies

Ensure these packages are installed:

```r
install.packages(c("ggplot2", "reshape2", "dplyr", "ComplexUpset", "stringr", "ggrepel"))
```

### Example Usage

```r
# Load trained XGBoost model
result <- xgboost_train(
  omics = c("Proteomics","Metabolomics","Microbiome", "Host"),
  feature = feature_df,
  seed = 123,
  top = 15,
  df.matrix = df_matrix,
  panel = "PreACLF",
  nfolds = 10,
  fecal.features = fecal_features,
  plasma.features = plasma_features,
  metadata = metadata_df
)

# Generate performance plots for top 10 features
Result.plots <- xgboost_eval(xgboost.object = result, top = 10)

# Access AUC plot
result$AUC.top.feature.curve
```

## Selection of the optimal feature-reduced model with  'Xgboost_best_model.GALAXY.R'

This workflow evaluates XGBoost feature-reduced models across the four clinical panels related to liver disease:**steatosis**, **inflammation**, **moderate fibrosis** and  **advanced fibrosis**. The objective is to identify the optimal number of features, selecting the model that used the fewest number of features while preserving at least 99% of the maximum average AUC.


##  Inputs

| Object Name            | Description                              | Notes                       |
|------------------------|------------------------------------------|-----------------------------|
| `xgboost.Fibrosis1`    | XGBoost model object for fibrosis ≥2. XGboost_train() output| Input to `xgboost_eval()`   |
| `xgboost.Fibrosis2`    | XGBoost model object for fibrosis ≥3. XGboost_train() output   | Input to `xgboost_eval()`   |
| `xgboost.Inflammation` | XGBoost model object for inflammation ≥2. XGboost_train() output | Input to `xgboost_eval()`   |
| `xgboost.steatosis`    | XGBoost model object for steatosis ≥2. XGboost_train() output    | Input to `xgboost_eval()`   |
| `xgboost_eval()`       | Custom function to evaluate model AUCs   | Sources(`xgboost_eval.GALAXY.R`) |



##  Outputs

| Object Name            | Description                                                      |
|------------------------|------------------------------------------------------------------|
| `best.auc.list`        | Nested list of models meeting ≥99% AUC and minimum feature number|
| `best.auc.df`          | Final dataframe with most efficient models per condition         |


##  Final Output Format (`best.auc.df`)

| Column        | Description                                 |
|---------------|---------------------------------------------|
| `omic`        | Omic data type (e.g., transcriptomics)       |
| `AUC.cv`      | Cross-validated AUC                         |
| `sd`          | Standard deviation of AUC                   |
| `panel`       | Clinical condition panel  (i.e.S>=1,I>=2,F>=2, F>=3 )|
| `no.features` | Number of features used                     |


## Evaluating optimal features-reduced model on the holdout test set with 'XGboost_validation()'

## Overview
This workflow validates the previously selected  feature-reduced  models for **steatosis**, **inflammation**, **moderate fibrosis** and  **advanced fibrosis** on the holdout test set. It evaluates the models on test data and calculates performance metrics such as AUC.

## Workflow
1. **Load Required Functions**  
  - Sources `xgboost_eval.GALAXY.R`, which contains the `xgboost_eval()` function.

2. **Load Data**  
- Holdout test set (`df.matrix.test`) and metadata test (`metadata`).
- Loads pre-trained XGBoost models for fibrosis, inflammation, and steatosis (output of xgboost_train() function).
- load output of XGboost_best_model.GALAXY.R (best.auc.df)

3. **Feature Selection for Testing**  
- Extracts relevant features (`fibrosis`, `kleiner_numeric`, `inflam_numeric`, `steatosis_numeric`) from metadata.
- Recategorizes variables based on thresholds for fibrosis, inflammation, and steatosis.

4. **Validation Function (`test.validation()`)**  
- Extracts the best models' features.
   - Evaluates models on test data and computes AUC.
   - Returns model validation metrics.

5. **Model Validation Execution**  
   - Calls `test.validation()` for each liver disease histopathologycal state (e.g.fibrosis (F>=2)) on the test data.

## Inputs and Outputs

### Inputs
| Variable | Type  | Description |
|----------|------------|------|
| `xgboost.model` | Model | XGBoost model object. XGboost_train() output |
| `df.matrix.test` | DataFrame     | Data matrix with omics features (columns) and samples (rows). **Note:** each feature must be prefixed with its corresponding omic type (e.g.,"Proteomics.P10643") 
| `metadata` | Dataframe  | test set metadata containing sample information  |
| `bestmodels.filt.AUC` | Dataframe |  Pre-filtered list of best models based on AUC. XGboost_best_model.GALAXY.R output |

### Outputs
| Variable | Type |  Description |
|----------|------------|------|
| `df.prediction.list` | List |  Predictions on the  test samples |
| `bestmodels.filt` | Dataframe| Filtered best models with validation AUC |


### Dependencies

- Ensure to source `xgboost_eval.GALAXY.R` for model evaluation.
- Ensure these packages are installed:

```r
install.packages(c( "Metrics", "dplyr"))
```


## Usage
1. Ensure the `xgboost_eval.GALAXY.R` script is available.
2. Load the required datasets.
3. Run `test.validation()` with appropriate parameters.

## Output
- A structured validation of the best models based on AUC.
- Final dataset (`bestmodels.filt`) containing optimal models for fibrosis, inflammation, and steatosis and its performance on the holdout test set.


# 3. Cox regression analysis in the GALA-ALD

To investigates the association between omics features and clinical outcomes of ALD patients (decompensation, mortality, and infection), we perform Cox proportional hazards regression analysis using `cox_regression.R`.

### Dependencies

```r
install.packages(c("tidyverse", "data.table", "survival"))
```

### Inputs

| Folder | Description |
|--------|-------------|
| `data/omic_data/` | Data matrix (TSV format) for omics datasets (e.g., cytokines, SNPs, microbiome) |

Each file is expected to:
- Have samples in rows (sample IDs as rownames)
- Have omic features in columns

### ️Analysis Overview

- Cox regression is run for **each feature** in each omics dataset.
- Each model is adjusted using clinical confounders from `GALA_ALD.confounder.rds`.
- The outcomes (decompensation, all-cause mortality, and infection) are tested separately.

### Outputs

| Folder | File Format | Description |
|--------|-------------|-------------|
| `out/cox_regression/` | `.tsv` | One file per input omics dataset per outcome tested |

Each file contains:
| Column | Description |
|--------|-------------|
| `name` | Feature name |
| `hz`   | Estimated hazard ratio |
| `p`    | p-value of the feature in the Cox model |
| `outcome` | Outcome type (decompensation, mortality, infection) |

# 4. Prognostic model construction for clinical outcomes in the GALA-ALD cohort

To construct prediction models for the clinical outcomes in the GALA-ALD, we employ randomforest model in the `mlr` package and make models using `construct_prognostic_model.R`.

### Dependencies
```r
install.packages(c("tidyverse", "mlr", "survival", "glmnet", "rsample", "data.table"))
```

### Inputs

| Folder | Description |
|--------|-------------|
| `data/omic_data/` | Data matrix (tsv format) for omics datasets (e.g., cytokines, SNPs, microbiome) |

Each file is expected to:
- Have samples in rows (sample IDs as rownames)
- Have features in columns.


### ️ Analysis Workflow
**Model Construction**  
- Repeated 5x5-fold cross-validation using `mlr::surv.ranger`
- Feature selection based on permutation importance
- C-index calculated for each model

### Outputs

All output is saved in: `out/prognostic_model/`

Each file includes:
| List Index | Content |
|------------|---------|
| `predicted_score` | Mean predicted risk score per sample across CV |
| `cindex` | C-index using the mean prediction |
| `model` | Full trained model (`mlr` object) |
| `surv` | Survival object used (`Surv`) |
| `d` | Original feature matrix (`d`) |
| `d2` | Modified matrix with `day` and `cond` (`d2`) |
| `cindex_fold` | C-index per fold and feature number (data frame) |
| `feature_rank` | Feature rankings per CV repetition and fold |
| `res` | Prediction results per fold and feature subset |

# 5. Summarize feature selection results for prognostic models

This script, `check_feature_selection_res.R`, summarizes the results of feature selection from prognostic survival models constructed using `construct_prognostic_model.R`. It identifies the smallest number of features that retain at least 99% of the best model performance (based on mean C-index), and extracts the corresponding feature names and performance metrics.

### Dependencies
```r
install.packages(c("tidyverse"))
```

### Inputs

| Location | Description |
|----------|-------------|
| `out/prognostic_model/*.rds` | Output from `construct_prognostic_model.R` containing full model results for each omics and outcome |

Each model `.rds` file is expected to contain:
- `cindex_fold`: C-index results for different feature numbers across CV folds
- `feature_rank`: Selected feature lists per fold

### Outputs

All outputs are saved to: `out/rds/`

| File | Description |
|------|-------------|
| `selected_feature.rds` | Data frame listing selected feature names for each omics-outcome combination, along with best fold identifiers (`i`, `j`) |
| `df.selected_feature.rds` | Summary table with best number of features and corresponding mean C-index |


# 6. Evaluation of prognostic models in the GALA-ALD

This script evaluates the performance of prognostic models trained using omics data from the GALA-ALD cohort. It compares both:
- Full models using **all features**
- Reduced models using **selected top features** (based on cross-validated C-index)

Evaluation metrics include:
- Concordance Index (C-index)
- Time-dependent AUC (1, 3, 5 years)
- Net Reclassification Index (NRI)

### Dependencies
```r
install.packages(c("tidyverse", "survival", "glmnet"))
# Use Bioconductor or remotes for:
# timeROC: https://cran.r-project.org/package=timeROC
# nricens: https://cran.r-project.org/package=nricens
```

### Inputs
### 1. Model Results
| Location | Description |
|----------|-------------|
| `out/prognostic_model/*.rds` | Full model objects from `construct_prognostic_model.R` |

Each file contains:
- Predicted survival scores (`res`)
- Fold assignments (`rep`, `fold`)
- Feature set used (`fs = "all"` or selected number)

### 2. Metadata
| File | Description |
|------|-------------|
| `data/metadata/GALA_ALD.marker.rds` | Metadata table including `sampleID` used to link predictions to survival outcome |

### 3. Selected Feature Summary
| File | Description |
|------|-------------|
| `out/rds/selected_feature.rds` | Identifies best feature set for each omics-outcome pair (used for reduced model evaluation) |

### Outputs

All results are saved in: `out/`

| File | Description |
|------|-------------|
| `evaluation_metrics.full_model_with_each_repetition.rds / .tsv` | Metrics for full models using all features |
| `evaluation_metrics.reduced_model_with_each_repetition.rds / .tsv` | Metrics for reduced models using selected features |

Each file contains:

| Column | Description |
|--------|-------------|
| `omics` | Data type (e.g., cytokines, SNPs) |
| `outcome` | Clinical outcome (decompensation, all-cause mortality, and infection) |
| `rep` | Repetition index (1 to 5) |
| `Cindex`, `AUC.xyr`, `NRI.xyr`, `n_total` | Evaluation metrics |
| `best_num` | (Reduced model only) Number of selected features |
