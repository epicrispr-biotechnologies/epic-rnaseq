process DESEQ2_REPORT {
    label "process_medium"

    // (Bio)conda packages have intentionally not been pinned to a specific version
    // This was to avoid the pipeline failing due to package conflicts whilst creating the environment when using -profile conda
    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/mulled-v2-8849acf39a43cdd6c839a369a74c0adc823e2f91:ab110436faf952a33575c64dd74615a84011450b-0' :
        '075615082992.dkr.ecr.us-west-2.amazonaws.com/deseq2_report:latest' }"

    input:
    val star_salmon_files
    val tx2gene
    path samplesheet
    path star_salmon_dir

    output:
    path "ls.txt", emit: ls_file
    path "deseq2_report.html", emit: report
    path "*.tsv", emit: deseq2_tsv

    script:
    def mypwd = new File('.').absolutePath
    """
    ls -l ${star_salmon_dir} > "ls.txt"
    echo "${moduleDir}"
    echo \$PATH
    which gtf2bed
    ls //nextflow-bin/
    head ${samplesheet}
    realpath ${samplesheet}
    echo \$PWD
    echo "'\${PWD}'"
    echo "'mypwd${mypwd}'"
    echo "Staging path: ${task.workDir}"
    Rscript -e "rmarkdown::render('//nextflow-bin/deseq2.rmd', output_file = 'deseq2_report.html', output_dir = '\${PWD}', intermediates_dir = '\${PWD}', params = list(study = 'RNAseq_1M', metadata_path = '\${PWD}/${samplesheet}', salmon_path = '\${PWD}/${star_salmon_dir}', deg_path = '\${PWD}/RNAseq_1M_degs.tsv'))"

    """
}
