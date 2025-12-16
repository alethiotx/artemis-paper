/**
 * Extract molecule data from ChEMBL database
 *
 * Queries the ChEMBL database for molecule information and exports deduplicated
 * records to CSV. ChEMBL is a large-scale bioactivity database containing
 * drug-like molecules with bioactivity data. The process uses pystow for
 * caching downloaded data.
 *
 * Note: Requires 200GB disk space for database download and processing.
 * Cache is stored in work directory (S3-backed via Fusion).
 */
process chembl {
  disk 200.GB
  label 'chembl'
  
  publishDir params.outdir + '/data/chembl', mode: 'copy', overwrite: false
  
  output:
    path "${params.chembl_version}/molecules.csv"
  
  script:
  """
  # Set ChEMBL cache to work directory (S3-backed via Fusion)
  export PYSTOW_HOME=\${PWD}/pystow_cache
  
  mkdir -p ${params.chembl_version}
  chembl.py "${params.chembl_version}"
  """
}

/**
 * Download and process MeSH tree structure
 *
 * Retrieves MeSH (Medical Subject Headings) hierarchical tree structure from
 * the National Library of Medicine (NLM) and serializes it to a pickle file.
 * MeSH provides a controlled vocabulary for indexing biomedical literature
 * and is used for disease classification and terminology standardization.
 *
 * The tree structure represents hierarchical relationships between medical
 * terms, enabling navigation from broad to specific concepts.
 */
process mesh {
  label 'mesh'
  
  publishDir params.outdir + '/data/mesh', mode: 'copy', overwrite: false
  
  output:
    path "*.pkl"
  
  script:
  """
  mesh.py "${params.mesh_file_base}"
  """
}