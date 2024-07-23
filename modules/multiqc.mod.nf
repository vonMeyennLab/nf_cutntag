#!/usr/bin/env nextflow
nextflow.enable.dsl=2


/* ========================================================================================
    PROCESSES
======================================================================================== */
process MULTIQC {

	label 'multiQC'

	container 'docker://josousa/multiqc:1.22.3'

	input:
		path(file)
		val(outputdir)
		val(multiqc_args)

	output:
		path "*html", emit: html

		publishDir "$outputdir/qc", mode: "link", overwrite: true

	script:
		"""
		multiqc ${multiqc_args} --filename multiqc_report.html .
		"""
}
