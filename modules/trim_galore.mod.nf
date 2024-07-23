#!/usr/bin/env nextflow
nextflow.enable.dsl=2


/* ========================================================================================
    PROCESSES
======================================================================================== */
process TRIM_GALORE {

	label 'trimGalore'
	tag "$name" // Adds name to job submission

	container 'docker://josousa/trimgalore:trimgalore-0.6.10-cutadapt-4.9'

	input:
		tuple val(name), path(reads)
		val(outputdir)
		val(trim_galore_args)

	output:
		tuple val(name), path("*fq.gz"), 			 emit: reads
		path "*trimming_report.txt", optional: true, emit: report
			
		publishDir "$outputdir/unaligned/fastq", mode: "link", overwrite: true, pattern: "*fq.gz"
		publishDir "$outputdir/unaligned/logs",  mode: "link", overwrite: true, pattern: "*trimming_report.txt"

	script:
	
		/* ==========
			Paired-end
		========== */
		if (reads instanceof List){
			trim_galore_args += " --paired "
		}

		"""
		trim_galore --cores ${task.cpus} ${trim_galore_args} ${reads}
		"""
}
