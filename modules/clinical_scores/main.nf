/**
 * Compute clinical scores for drug targets across disease indications
 *
 * Queries ChEMBL for drug-target-disease associations and clinical trial data
 * to generate scored target lists for seven disease indications. Combines
 * clinical trial phase information with trial recency to prioritize the most
 * relevant therapeutic targets.
 *
 * Computes scores for:
 * - Breast cancer (Breast Neoplasms)
 * - Lung cancer (Lung Neoplasms)
 * - Prostate cancer (Prostatic Neoplasms)
 * - Melanoma (Skin Neoplasms)
 * - Bowel cancer (Intestinal Neoplasms)
 * - Type 2 diabetes (Diabetes Mellitus, Type 2)
 * - Cardiovascular disease (Cardiovascular Diseases)
 *
 * Outputs are dated directories containing CSV files with target genes,
 * phase scores, and clinical scores for each indication.
 */
process clinical_scores {
  label 'clinical_scores'
  
  publishDir params.outdir + '/data', mode: 'copy', overwrite: false
  
  output:
    path '*', emit: csv
  
  script:
  """
  mkdir -p clinical_scores
  clinical_scores.py
  """
}