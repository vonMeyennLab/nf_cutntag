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


process BEDTOOLS_GENOMECOV_NORM {

    label 'bedtools'
    tag "$bam" // Adds name to job submission

    container 'docker://staphb/bedtools:2.31.1'

    input:
        path bam
		path metrics
        val outputdir
        val bedtools_genomecov_args

    output:
        path "*bedgraph", emit: bedgraph
        publishDir "$outputdir/aligned/bedgraph/normalized/", mode: "link", overwrite: true

    script:
        """

		# get number of mapped reads from Picard MarkDupl metrics file
		mapped=\$(grep -A2 " METRICS CLASS" ${metrics} | grep -v "#" | cut -f 3 | grep -v READ_PAIRS_EXAMINED)
		# get norm factor, by dividing the number of mapped reads by 1M
		mapped=\$(awk -v var=\$mapped 'BEGIN { printf "%.2f", var/1000000 }')

		# convert bam to bedgraph
        bedtools genomecov -scale \${mapped} ${bedtools_genomecov_args} -ibam ${bam} > ${bam}.normalized.bedgraph
        
        for file in *.bam.normalized.bedgraph; do
            mv "\$file" "\${file/.bam.normalized.bedgraph/.normalized.bedgraph}"
        done
        """
}
