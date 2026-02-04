#!/usr/bin/env nextflow
nextflow.enable.dsl=2


/* ========================================================================================
    DEFAULT PARAMETERS
======================================================================================== */
params.bam_output = true // Setting if the bam file should be published


/* ========================================================================================
    PROCESSES
======================================================================================== */
process BOWTIE2 {

		label 'bowtie2'
		tag "$name" // Adds name to job submission

		container 'docker://josousa/bowtie2:2.5.4'

	input:
		tuple val(name), path(reads)
		val(outputdir)
		val(bowtie2_args)

	output:
		tuple path "*bam", path "*stats.txt"
		//path "*bam",  	   emit: bam
		//path "*stats.txt", emit: stats

		publishDir "$outputdir/aligned/bam",  mode: "link", overwrite: true, pattern: "*bam", enabled: params.bam_output
		publishDir "$outputdir/aligned/logs", mode: "link", overwrite: true, pattern: "*stats.txt"

	script:
		/* ==========
			File names
		========== */
		readString = ""
		if (reads instanceof List) {
			readString = "-1 " + reads[0] + " -2 " + reads[1]
			bowtie2_args += " --no-discordant --no-mixed "
		}
		else {
			readString = "-U " + reads
		}

		/* ==========
			Index
		========== */
		index = params.genome["bowtie2"]

		/* ==========
			Basename
		========== */
		bowtie_name = name + "_" + params.genome["name"]

		"""
		bowtie2 -x ${index} -p ${task.cpus} ${bowtie2_args} ${readString} 2>${bowtie_name}_bt2_stats.txt | samtools view -@ ${task.cpus} -bS -F 4 -F 8 -F 256 -> ${bowtie_name}_bt2.bam
		"""
}
