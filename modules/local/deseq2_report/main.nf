process DESEQ2_REPORT {
    label "process_medium"

    // (Bio)conda packages have intentionally not been pinned to a specific version
    // This was to avoid the pipeline failing due to package conflicts whilst creating the environment when using -profile conda
    conda "${moduleDir}/environment.yml"
    container "'075615082992.dkr.ecr.us-west-2.amazonaws.com/deseq2_report:latest"

    input:
    val star_salmon_files
    val tx2gene       
    path samplesheet
    path star_salmon_dir
    val user
    val study

    output:
    path "*deseq2_report.html", emit: report
    path "*.tsv", emit: deseq2_tsv

    script:
    """
    Rscript -e "rmarkdown::render(
        '//nextflow-bin/deseq2.rmd',
        output_file = '${study}_deseq2_report.html',
        output_dir = '\${PWD}',
        intermediates_dir = '\${PWD}',
        params = list(
            study = '${study}',
            metadata_path = '\${PWD}/${samplesheet}',
            salmon_path = '\${PWD}/${star_salmon_dir}',
            deg_path = '\${PWD}/${study}_degs.tsv',
            user = '${user}'
        )
    )"

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        r-base: \$(echo \$(R --version 2>&1) | sed 's/^.*R version //; s/ .*\$//')
        bioconductor-deseq2: \$(Rscript -e "library(DESeq2); cat(as.character(packageVersion('DESeq2')))")
    END_VERSIONS
    """

    stub:
    """
    touch deseq2_report.html
    touch degs.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        r-base: \$(echo \$(R --version 2>&1) | sed 's/^.*R version //; s/ .*\$//')
        bioconductor-deseq2: \$(Rscript -e "library(DESeq2); cat(as.character(packageVersion('DESeq2')))")
    END_VERSIONS
    """
}
