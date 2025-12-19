
/**
 * Compute drug target predictions for a specific parameter combination
 *
 * Runs the compute.py script to generate predictions across multiple disease
 * indications (breast, lung, bowel, prostate, melanoma, diabetes, cardiovascular).
 * Outputs include target probabilities, SABCS overlaps, and training data.
 *
 * @param kg Knowledge graph name (hetionet_full, biokg, openbiolink, primekg)
 * @param ct_unique Clinical target filter threshold
 * @param probs Probability threshold for predictions
 * @param p_genes Pathway gene inclusion parameter
 * @param iteration Cross-validation iteration number (1-10)
 */
process compute {
  tag "${kg}-${ct_unique}-${probs}-${p_genes}-${iteration}"
  
  // publishDir params.outdir + '/figs/predictions', mode: 'copy', pattern: 'targets/*.csv'
  // publishDir params.outdir + '/figs/predictions', mode: 'copy', pattern: 'training/*.csv'
  
  input:
    tuple val(kg), val(ct_unique), val(probs), val(p_genes), val(iteration)
  
  output:
    path('predictions/*.csv'), emit: predictions
    path('sabcs/*.csv'), emit: sabcs
    path('targets/*.csv'), emit: targets
    path('training/*.csv'), emit: training
  
  script:
  """
  mkdir -p training targets predictions sabcs
  
  compute.py \\
    ${kg} \\
    ${ct_unique} \\
    ${probs} \\
    ${p_genes} \\
    ${iteration} \\
    ${params.scores_date}
  """
}

/**
 * Combine and visualize prediction results across all parameter combinations
 *
 * Aggregates prediction outputs from all compute processes, calculates overlap
 * metrics (sensitivity, specificity), and generates visualization heatmaps
 * grouped by knowledge graph and disease indication.
 *
 * @param csv Collection of prediction CSV files from all compute iterations
 */
process pr_combine {
  publishDir params.outdir + '/figs/predictions', mode: 'copy'
  
  input:
    path(csv)
  
  output:
    path 'plots/*'
    path 'data/*.pickle'
    path 'data/all.csv', emit: combined_csv
  
  script:
  """
  mkdir -p plots/indications plots/kgs plots/heatmaps data
  
  combine.py ${csv}
  """
}

/**
 * Aggregate target probabilities across all indications and parameter combinations
 *
 * Processes all prediction CSV files for all disease indications, organizing
 * target probabilities by indication, knowledge graph, clinical threshold, probability
 * threshold, and pathway gene parameters. Outputs a single serialized pickle file
 * containing all indications for downstream baseline analysis.
 *
 * @param csv Collection of all target prediction CSV files
 */
process targets {
  label 'targets'
  publishDir params.outdir + '/figs/predicted_targets', mode: 'copy'
  
  input:
    path(csv)
  
  output:
    path('all_targets.pickle'), emit: pickle
  
  script:
  """
  targets.py ${csv}
  """
}

/**
 * Aggregate training set labels across all indications and parameter combinations
 *
 * Processes all training set CSV files for all disease indications, organizing
 * binary labels by indication, knowledge graph, clinical threshold, probability
 * threshold, and pathway gene parameters. Outputs a single serialized pickle file
 * containing all training labels for downstream analysis.
 *
 * @param csv Collection of all training set CSV files
 */
process training_sets {
  label 'training'
  publishDir params.outdir + '/figs/training_sets', mode: 'copy'
  
  input:
    path(csv)
  
  output:
    path('all_training_sets.pickle'), emit: pickle
  
  script:
  """
  training_sets.py ${csv}
  """
}

/**
 * Compute baseline prediction statistics and visualize parameter effects
 *
 * Analyzes aggregated target probabilities to compute baseline statistics
 * (percentage of targets above Random Forest thresholds) across different
 * parameter configurations for each indication. Generates heatmaps showing
 * the effect of clinical target filtering and probability thresholds.
 *
 * @param pickle Serialized target probability dictionary for all indications
 */
process baselines {
  label 'baselines'
  publishDir params.outdir + '/figs/baselines', mode: 'copy'
  
  input:
    path(pickle)
    path(combined_csv)
  
  output:
    path 'plots/*'
    path 'data/*'
  
  script:
  """
  mkdir -p plots data
  
  baselines.py ${pickle} ${combined_csv}
  """
}

/**
 * Analyze overlap between predictions and SABCS breast cancer targets
 *
 * Compares predicted drug targets against the SABCS (San Antonio Breast Cancer
 * Symposium) dataset to evaluate prediction accuracy for known breast cancer
 * targets. Generates heatmaps showing overlap percentages across different
 * knowledge graphs and Random Forest probability thresholds.
 *
 * @param csv Collection of SABCS overlap CSV files from compute processes
 */
process sabcs {
  publishDir params.outdir + '/figs/sabcs', mode: 'copy'
  
  input:
    path(csv)
  
  output:
    path 'plots/*'
    path 'data/*'
  
  script:
  """
  mkdir -p plots data
  
  sabcs.py ${csv}
  """
}

/**
 * Generate consensus predictions across knowledge graphs for SABCS targets
 *
 * Computes consensus predictions by averaging target probabilities across all
 * knowledge graphs (Hetionet, BioKG, OpenBioLink, PrimeKG) and iterations.
 * Creates a clustered heatmap showing the agreement between different KGs
 * for SABCS breast cancer targets.
 *
 * @param pickle Target probability dictionary for all indications
 */
process sabcs_consensus {
  label 'sabcs_consensus'
  publishDir params.outdir + '/figs/sabcs_consensus', mode: 'copy'
  
  input:
    path(pickle)
  
  output:
    path 'plots/*'
    path 'data/*'
  
  script:
  """
  mkdir -p plots data
  
  consensus_sabcs.py ${params.scores_date} ${pickle}
  """
}
