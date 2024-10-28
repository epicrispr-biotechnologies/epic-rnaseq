process DESEQ2_REPORT {
    label "process_medium"

    // (Bio)conda packages have intentionally not been pinned to a specific version
    // This was to avoid the pipeline failing due to package conflicts whilst creating the environment when using -profile conda
    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/mulled-v2-8849acf39a43cdd6c839a369a74c0adc823e2f91:ab110436faf952a33575c64dd74615a84011450b-0' :
        'biocontainers/mulled-v2-8849acf39a43cdd6c839a369a74c0adc823e2f91:ab110436faf952a33575c64dd74615a84011450b-0' }"

    input:
    val star_salmon_files
    path tx2gene
    path samplesheet
    path star_salmon_dir

    output:
    path "ls.txt", emit: ls_file
    path "deseq2.html", emit: report
    path "*.tsv", emit: deseq2_tsv

    script:
    """
    ls -l ${star_salmon_dir} > "ls.txt"
    Rscript -e "rmarkdown::render('bin/deseq2.rmd', params = list(study = 'RNAseq_1M', metadata_path = ${samplesheet}, salmon_path = ${star_salmon_dir}, deg_path = 'RNAseq_1M_degs.tsv'))"

    """
}
