#!/usr/bin/env nextflow
nextflow.enable.dsl=2


/* ========================================================================================
    PROCESSES
======================================================================================== */
process MARK_DUPLICATES {	

	label 'picard'
	tag "$bam" // Adds name to job submission

	container 'docker://agrf/picard:2.18.27'

	input:
		tuple val(name), path(bam)
		//path(bam)
		val(outputdir)
		val(mark_duplicates_args)

	output:
		tuple val(name), path("*bam"), emit: bam
		tuple val(name), path("*txt"), emit: metrics

		publishDir "$outputdir/aligned/bam",              mode: "link", overwrite: true, pattern: "*dupflag.bam"
		publishDir "$outputdir/aligned/bam/deduplicated", mode: "link", overwrite: true, pattern: "*dedup.bam", enabled: params.bam_output
		publishDir "$outputdir/aligned/logs",             mode: "link", overwrite: true, pattern: "*txt"

	script:
		/* ==========
			Basename
		========== */
		base_name = bam.toString() - ".sorted.bam"

		/* ==========
			Mark duplicates
		========== */
		if ((mark_duplicates_args =~ /.*REMOVE_DUPLICATES=true.*/)) {
			output_name = base_name + ".sorted.dedup.bam"
		} else if (!(mark_duplicates_args =~ /.*REMOVE_DUPLICATES=true.*/) && (mark_duplicates_args =~ /.*REMOVE_SEQUENCING_DUPLICATES=true.*/)) {
			output_name = base_name + ".sorted.seqdedup.bam"
		} else {
			output_name = base_name + ".sorted.dupflag.bam"
		}

		"""
		picard MarkDuplicates INPUT=${bam} OUTPUT=${output_name} ASSUME_SORTED=true METRICS_FILE=${base_name}.MarkDuplicates.metrics.txt ${mark_duplicates_args}
		"""
}
