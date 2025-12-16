/**
 * Perform cross-validation for a specific knowledge graph, indication, and use case
 *
 * Executes the appropriate CV script (binary, multiclass, or regression) based on
 * the use case parameter. Each script performs k-fold cross-validation with multiple
 * classifiers/regressors and outputs performance metrics.
 *
 * Binary: AUROC, accuracy, sensitivity, specificity with RandomForest and SVM
 * Multiclass: AUROC, accuracy across 2, 3, and 5 bins with RandomForest and SVM
 * Regression: R² with LinearRegression and KNeighborsRegressor
 *
 * @param kg Knowledge graph name (hetionet_full, biokg, openbiolink, primekg)
 * @param indication Disease indication (breast, lung, bowel, prostate, melanoma, diabetes, cardiovascular)
 * @param usecase Analysis type (binary, multiclass, regression)
 */
process cv {
  tag "${kg}-${indication}-${usecase}"
  label 'cv'

  publishDir params.outdir + '/figs/cv', mode: 'copy'
  
  input:
    tuple val(kg), val(indication), val(usecase)
  
  output:
    path 'data/*.csv', emit: csv
  
  script:
  """
  if [[ "${usecase}" == "binary" ]]; then
    binary.py ${kg} ${indication} ${params.scores_date}
  elif [[ "${usecase}" == "multiclass" ]]; then
    multiclass.py ${kg} ${indication} ${params.scores_date}
  elif [[ "${usecase}" == "regression" ]]; then
    regression.py ${kg} ${indication} ${params.scores_date}
  else
    echo "Unknown usecase parameter: ${usecase}"
    exit 1
  fi
  """
}

/**
 * Combine and visualize cross-validation results across all configurations
 *
 * Aggregates CV results from all knowledge graphs, indications, and use cases.
 * Generates standardized plots comparing performance metrics across different
 * configurations. Creates heatmaps and comparison plots grouped by use case,
 * knowledge graph, and disease indication.
 *
 * @param csv Collection of CV result CSV files from all cv processes
 */
process cv_combine {
  publishDir params.outdir + '/figs/cv', mode: 'copy'
  
  input:
    file csv
  
  output:
    path('data/*.csv')
    path('plots/*.png')
  
  script:
  """
  combine.py ${csv}
  """
}

/**
 * Analyze cross-validation performance across feature subsampling ranges
 *
 * Performs feature subsampling analysis to understand how the number of features
 * affects model performance. Uses a fixed configuration (hetionet_full, breast)
 * with multiple iterations to evaluate performance at different feature counts.
 * This helps determine the optimal number of features needed for predictions.
 *
 * @param subsample Feature subsample size to evaluate
 */
process range {
  tag "${subsample}"
  label 'cv_range'

  publishDir params.outdir + '/figs/cv_range', mode: 'copy'
  
  input:
    val subsample
  
  output:
    path 'data/*.csv', emit: csv
  
  script:
  """
  range.py ${subsample} ${params.scores_date}
  """
}

/**
 * Combine and visualize feature range analysis results
 *
 * Aggregates results from feature subsampling experiments across all subsample
 * sizes. Computes median trends with polynomial smoothing (degree 7) and generates
 * plots showing how AUROC performance changes with the number of features.
 * Helps identify the point of diminishing returns for feature inclusion.
 *
 * @param csv Collection of feature range CSV files from all range processes
 */
process range_combine {
  publishDir params.outdir + '/figs/cv_range', mode: 'copy'
  
  input:
    file csv
  
  output:
    path('data/*.csv')
    path('plots/*.png')
  
  script:
  """
  range_combine.py ${csv}
  """
}
