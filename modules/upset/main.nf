/**
 * Generate UpSet plots for clinical targets and pathway genes
 *
 * Creates UpSet plots showing overlaps and unique sets across disease indications
 * for both clinical targets and pathway genes. Generates five visualizations:
 * - All clinical targets overlap
 * - Approved clinical targets only
 * - Unique clinical targets per indication
 * - Top 100 pathway genes overlap
 * - Top 300 pathway genes overlap
 *
 * UpSet plots are useful for visualizing set intersections and understanding
 * which targets or genes are shared across multiple disease indications versus
 * those that are indication-specific.
 */
process upset {
  publishDir params.outdir + '/figs/upset', mode: 'copy'
  
  output:
    path '*.png'
  
  script:
  """
  upset.py ${params.scores_date}
  """
}
