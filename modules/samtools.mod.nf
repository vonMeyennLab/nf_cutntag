#!/usr/bin/env nextflow
nextflow.enable.dsl=2


/* ========================================================================================
    DEFAULT PARAMETERS
======================================================================================== */
params.bam_output = true // Setting if the bam file should be published


/* ========================================================================================
    PROCESSES
======================================================================================== */

// SAMTOOLS_SORT
process SAMTOOLS_SORT {	
    
	label 'samtools'
	tag "$bam" // Adds name to job submission

	container 'docker://staphb/samtools:1.20'

	input:
		path(bam)
		val(outputdir)
		val(samtools_sort_args)

	output:
		path "*bam", emit: bam
		publishDir "$outputdir/aligned/bam", mode: "link", overwrite: true, enabled: params.bam_output

    script:
		/* ==========
			Basename
		========== */
		basename = bam.toString() - ".bam"

		"""
		samtools sort --threads ${task.cpus-1} ${samtools_sort_args} ${bam} -o ${basename}.sorted.bam
		"""
}

// SAMTOOLS_INDEX
process SAMTOOLS_INDEX {	
    
	label 'samtools'
	tag "$bam" // Adds name to job submission

	container 'docker://staphb/samtools:1.20'

	input:
		path(bam)
		val(outputdir)
		val(samtools_index_args)

	output:
		path "*.bai", emit: bai
		publishDir "$outputdir/aligned/bam", mode: "link", overwrite: true

    script:

		"""
		samtools index -@ ${task.cpus-1} ${samtools_index_args} ${bam}
		"""
}
