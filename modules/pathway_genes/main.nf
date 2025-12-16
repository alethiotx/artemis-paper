/**
 * Extract pathway genes for disease indications
 *
 * Queries the Reactome pathway database to identify genes involved in biological
 * pathways relevant to each disease indication. Pathway genes represent biological
 * mechanisms and processes associated with disease pathophysiology.
 *
 * Extracts pathway genes for seven indications:
 * - Breast Cancer
 * - Lung Cancer
 * - Bowel Cancer
 * - Prostate Cancer
 * - Melanoma
 * - Diabetes Mellitus Type 2
 * - Cardiovascular Disease
 *
 * Outputs are dated directories containing CSV files with gene lists for each
 * indication. These genes can be used as features or filters in prediction models.
 */
process pathway_genes {
  label 'pathway_genes'
  
  publishDir params.outdir + '/data', mode: 'copy', overwrite: false
  
  output:
    path '*', emit: csv
  
  script:
  """
  mkdir -p pathway_genes
  pathway_genes.py
  """
}