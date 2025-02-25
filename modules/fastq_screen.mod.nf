#!/usr/bin/env nextflow
nextflow.enable.dsl=2


/* ========================================================================================
    DEFAULT PARAMETERS
======================================================================================== */
params.fastq_screen_conf = '/cluster/work/nme/software/config/fastq_screen.conf'


/* ========================================================================================
    PROCESSES
======================================================================================== */
process FASTQ_SCREEN {

	label 'fastqScreen'
	tag "$name" // Adds name to job submission

	container 'docker://josousa/fastq_screen:0.15.3'
	
	input:
		tuple val(name), path(reads)
		val(outputdir)
		val(fastq_screen_args)

	output:
		path "*html", emit: html
		path "*txt",  emit: report
		
		publishDir "$outputdir/qc/fastq_screen", mode: "link", overwrite: true

	script:
		/* ==========
			Paired-end
		========== */
		if (reads instanceof List) {
			reads = reads[0]
		}

		"""
		fastq_screen --threads ${task.cpus} --conf ${params.fastq_screen_conf} ${fastq_screen_args} ${reads}
		"""
}
