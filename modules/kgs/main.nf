/**
 * Generate knowledge graph overview and feature analysis
 *
 * Executes a Jupyter notebook that analyzes and compares different knowledge graphs
 * used in the pipeline (Hetionet, BioKG, OpenBioLink, PrimeKG). The notebook
 * generates statistics about graph structure, node/edge types, and feature
 * distributions. Outputs an HTML report with interactive visualizations and
 * summary tables.
 *
 * The notebook is executed and converted to HTML format using nbconvert with
 * the lab template for enhanced visualization support.
 *
 * @param ipynb Jupyter notebook file containing KG analysis code
 */
process kgs_overview {
  label 'kgs'
  
  publishDir params.outdir + '/data/kgs/features', mode: 'copy'
  
  input:
    path ipynb
  
  output:
    path '*.html', emit: html
  
  script:
  """
  jupyter nbconvert \\
    --to html \\
    --execute \\
    --template lab \\
    --output kgs.html \\
    --output-dir . \\
    ${ipynb}
  """
}
