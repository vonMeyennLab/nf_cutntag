#!/usr/bin/env nextflow
nextflow.enable.dsl=2


/* ========================================================================================
	PROCESSES
======================================================================================== */
process BEDTOOLS_GENOMECOV {

    label 'bedtools'
    tag "$bam" // Adds name to job submission

    container 'docker://staphb/bedtools:2.31.1'

    input:
        path bam
        val outputdir
        val bedtools_genomecov_args

    output:
        path "*bedgraph", emit: bedgraph
        publishDir "$outputdir/aligned/bedgraph/raw/", mode: "link", overwrite: true

    script:
        """
        bedtools genomecov ${bedtools_genomecov_args} -ibam ${bam} > ${bam}.bedgraph
        
        for file in *.bam.bedgraph; do
            mv "\$file" "\${file/.bam.bedgraph/.bedgraph}"
        done
        """
}
