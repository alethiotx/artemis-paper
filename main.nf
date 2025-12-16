/**
 * Knowledge Graph-Based Drug Target Discovery Pipeline
 *
 * Main Nextflow workflow orchestrating knowledge graph analysis, cross-validation,
 * and drug target prediction across multiple disease indications.
 *
 * Pipeline Modes:
 * - kgs: Generate knowledge graph overview and feature analysis
 * - upset: Create UpSet plots for target/gene overlaps
 * - data: Download ChEMBL molecules and MeSH tree structure
 * - scores: Compute clinical scores and extract pathway genes
 * - cv: Perform cross-validation across KGs and indications
 * - predictions: Generate drug target predictions and consensus analysis
 *
 * Knowledge Graphs:
 * - Hetionet: Integrative biomedical knowledge graph
 * - BioKG: Biological knowledge graph
 * - OpenBioLink: Open biomedical link prediction dataset
 * - PrimeKG: Precision medicine knowledge graph
 *
 * Disease Indications:
 * - Breast, lung, bowel, prostate cancer
 * - Melanoma
 * - Type 2 diabetes
 * - Cardiovascular disease
 */

include { chembl; mesh } from './modules/data'
include { clinical_scores } from './modules/clinical_scores'
include { pathway_genes } from './modules/pathway_genes'
include { cv; cv_combine; range; range_combine } from './modules/cv'
include { compute; combine; targets; training_sets; baselines; sabcs; sabcs_consensus } from './modules/predictions'
include { upset } from './modules/upset'
include { kgs_overview } from './modules/kgs'

workflow {
  // ─── Knowledge Graphs ──────────────────────────────────────────────────────
  
  kgs = Channel.of(
    'hetionet',
    'biokg',
    'openbiolink',
    'primekg'
  )

  // ─── Knowledge Graph Overview ─────────────────────────────────────────────
  
  if (params.mode == 'kgs') {
    ipynb = Channel.fromPath('s3://alethiotx-artemis/data/kgs/features/kgs.ipynb')
    kgs_overview(ipynb)
  }

  // ─── UpSet Plot Generation ────────────────────────────────────────────────
  
  if (params.mode == 'upset') {
    upset()
  }

  // ─── External Data Download ───────────────────────────────────────────────
  
  if (params.mode == 'data') {
    chembl()
    mesh()
  }

  // ─── Clinical Scores & Pathway Genes ──────────────────────────────────────
  
  if (params.mode == 'scores') {
    clinical_scores()
    pathway_genes()
  }

  // ─── Cross-Validation ─────────────────────────────────────────────────────
  
  if (params.mode == 'cv') {
    // Disease indications
    indications = Channel.of(
      "breast",
      "lung",
      "bowel",
      "prostate",
      "melanoma",
      "diabetes",
      "cardiovascular"
    )
    
    // Analysis types: binary, multiclass (2/3/5 bins), regression
    usecases = Channel.of(
      'binary',
      'multiclass',
      'regression'
    )
    
    // Run CV across all KG x indication x use case combinations
    cv(
      kgs
        .combine(indications)
        .combine(usecases)
    )

    // Aggregate and visualize CV results
    cv_combine(
      cv.out
        .collect()
    )

    // Feature subsampling analysis (10,000 total features in Hetionet)
    subsample = Channel.of(30000, 20000, 10000, 5000, 1000, 500, 100, 50, 10, 5)
    
    range(
      subsample
    )

    // Combine feature range results with polynomial smoothing
    range_combine(
      range.out
        .collect()
    )
  }

  // ─── Drug Target Predictions ──────────────────────────────────────────────
  
  if (params.mode == 'predictions') {
    // Clinical trial filters: All, Unique per indication, Approved only
    ct_unique = Channel.of('All', 'Unique', 'Approved')
    
    // Random Forest probability thresholds
    probs = Channel.of(0.5, 0.6, 0.7, 0.8, 0.9)
    
    // Pathway gene feature counts: 0 (none), 100 (top), 300 (extended)
    p_genes = Channel.of(0, 100, 300)
    
    // Cross-validation iterations for stability
    iterations = Channel.of(1, 2, 3, 4, 5, 6, 7, 8, 9, 10)

    // Generate predictions across all parameter combinations
    compute(
      kgs
        .combine(ct_unique)
        .combine(probs)
        .combine(p_genes)
        .combine(iterations)
    )

    // Aggregate prediction overlaps and compute sensitivity/specificity
    combine(
      compute.out.predictions
        .collect()
    )

    // Aggregate target probabilities for all indications
    targets(
      compute.out.targets
        .collect()
    )

    // Aggregate training set labels for all indications
    training_sets(
      compute.out.training
        .collect()
    )

    // Analyze SABCS breast cancer target overlaps
    sabcs(
      compute.out.sabcs
        .collect()
    )

    // Compute baseline statistics by indication
    baselines(
      targets.out.pickle
    )

    // Generate consensus predictions across KGs for SABCS targets
    sabcs_consensus(
      targets.out.pickle
    )
  }

  // ─── Post Predictions calculations only (no need to rerun compute)──────────────────────────────────────────────

  if (params.mode == 'post_predictions') {

    // Compute baseline statistics by indication
    baselines(
      's3://alethiotx-artemis/figs/predicted_targets/all_targets.pickle'
    )

    // Generate consensus predictions across KGs for SABCS targets
    sabcs_consensus(
      's3://alethiotx-artemis/figs/predicted_targets/all_targets.pickle'
    )
  }

}