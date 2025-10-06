#!/usr/bin/env nextflow
nextflow.enable.dsl=2


/* ========================================================================================
    INPUT FILES
======================================================================================== */
params.input = null

if (!params.input) {
    error "Input not specified. Use --input to specify the input."
}

input_files = file(params.input)


/* ========================================================================================
    OUTPUT DIRECTORY
======================================================================================== */
params.outdir = false
if(params.outdir){
    outdir = params.outdir
} else {
    outdir = '.'
}


/* ========================================================================================
    SKIP STEPS
======================================================================================== */
params.skip_qc             = false
params.skip_fastq_screen   = true


/* ========================================================================================
    DEFAULT PARAMETERS
======================================================================================== */
params.fastq_screen_conf = '/cluster/work/nme/software/config/fastq_screen.conf' // FastQ Screen config file directory
params.genome            = 'Mus_musculus_GRCm39' // Default genome
params.single_end        = false    // Force to input files to be single-end


/* ========================================================================================
    PARAMETERS
======================================================================================== */
params.fastqc_args             = ''
params.fastq_screen_args       = ''
params.trim_galore_args        = ''
// Default Bowtie2 CutnTag parameters
params.bowtie2_args            = '--local --very-sensitive-local --minins 10 --maxins 700'
params.samtools_sort_args      = ''
params.samtools_index_args     = ''
params.mark_duplicates_args    = ''
params.bedtools_genomecov_args = ''
params.multiqc_args            = ''

fastqc_args             = params.fastqc_args
fastq_screen_args       = params.fastq_screen_args
trim_galore_args        = params.trim_galore_args
bowtie2_args            = params.bowtie2_args
samtools_sort_args      = params.samtools_sort_args
samtools_index_args     = params.samtools_index_args
mark_duplicates_args    = params.mark_duplicates_args
bedtools_genomecov_args = params.bedtools_genomecov_args
multiqc_args            = params.multiqc_args


/* ========================================================================================
    BOWTIE2 PARAMETERS
======================================================================================== */

// --no-unal
params.bowtie2_no_unal  = true
if (params.bowtie2_no_unal){
		bowtie2_args += " --no-unal "
	}
// suppress SAM records for unaligned reads


/* ========================================================================================
    PICARD MARKDUPLICATES PARAMETERS
======================================================================================== */

// --REMOVE_DUPLICATES
params.picard_markduplicates_remove_duplicates = false
if (params.picard_markduplicates_remove_duplicates){
        mark_duplicates_args += " REMOVE_DUPLICATES=true "
    }
// If true do not write duplicates to the output file instead of writing them with
// appropriate flags set.


// --REMOVE_SEQUENCING_DUPLICATES
params.picard_markduplicates_remove_sequencing_duplicates = false
if (params.picard_markduplicates_remove_sequencing_duplicates){
        mark_duplicates_args += " REMOVE_SEQUENCING_DUPLICATES=true "
    }
// If true remove 'optical' duplicates and other duplicates that appear to have arisen from
// the sequencing process instead of the library preparation process, even if
// REMOVE_DUPLICATES is false. If REMOVE_DUPLICATES is true, all duplicates are removed and
// this option is ignored.

// ---------------------------------------- //
// Condition to output the original sorted.bam file only if that file will be subsequently deduplicated

// With deduplication,    outputs:     *.sorted.bam and *.sorted.dedup.bam
// Without deduplication, outputs:     *.sorted.dupflag.bam
if ((params.mark_duplicates_args =~ /.*REMOVE_DUPLICATES=true.*/) || (params.mark_duplicates_args =~ /.*REMOVE_SEQUENCING_DUPLICATES=true.*/)){
        params.bam_output = true
    } else {
        params.bam_output = false
    }
// ---------------------------------------- //


/* ========================================================================================
    BEDTOOLS GENOMECOV PARAMETERS
======================================================================================== */

// -bga
params.bedtools_genomecov_bga = true
if (params.bedtools_genomecov_bga) {
	bedtools_genomecov_args += " -bga "
}
/* ==========
Report depth in BedGraph format, as '-bg'.
However with this option, regions with zero
coverage are also reported. This allows one to
quickly extract all regions of a genome with 0
coverage by applying: "grep -w 0$" to the output.
========== */


/* ========================================================================================
    GENOMES
======================================================================================== */
params.custom_genome_file = '' // Option to add a directory for a custom genome file

include { getGenome } from './modules/genomes.mod.nf' params(custom_genome_file: params.custom_genome_file)
genome = getGenome(params.genome)


/* ========================================================================================
    FILES CHANNEL
======================================================================================== */
include { makeFilesChannel; getFileBaseNames } from './modules/files.mod.nf'
file_ch = makeFilesChannel(input_files)


/* ========================================================================================
    WORKFLOW
======================================================================================== */
include { FASTQC }                  from './modules/fastqc.mod.nf'
include { FASTQC as FASTQC2 }       from './modules/fastqc.mod.nf'
include { FASTQ_SCREEN }            from './modules/fastq_screen.mod.nf'
include { TRIM_GALORE }             from './modules/trim_galore.mod.nf'
include { BOWTIE2 }                 from './modules/bowtie2.mod.nf'  params(bam_output: false, genome: genome)
include { SAMTOOLS_SORT }           from './modules/samtools.mod.nf' params(bam_output: params.bam_output)
include { SAMTOOLS_INDEX }          from './modules/samtools.mod.nf' params(bam_output: params.bam_output)
include { MARK_DUPLICATES }         from './modules/picard.mod.nf'   params(bam_output: params.bam_output)
include { BEDTOOLS_GENOMECOV_NORM } from './modules/bedtools.mod.nf'
include { BEDTOOLS_GENOMECOV }      from './modules/bedtools.mod.nf'
include { MULTIQC }                 from './modules/multiqc.mod.nf'

workflow {

    main:
        // QC conditional
        if (!params.skip_qc){ 
            FASTQC                          (file_ch, outdir, fastqc_args)
            
            // FastQ Screen conditional
            if (!params.skip_fastq_screen){ 
            FASTQ_SCREEN                    (file_ch, outdir, fastq_screen_args)
            }

            TRIM_GALORE                     (file_ch, outdir, trim_galore_args)
            FASTQC2                         (TRIM_GALORE.out.reads, outdir, fastqc_args)
            BOWTIE2                         (TRIM_GALORE.out.reads, outdir, bowtie2_args)
        } else {
            BOWTIE2                         (file_ch, outdir, bowtie2_args)
        }
        SAMTOOLS_SORT             (BOWTIE2.out.bam, outdir,samtools_sort_args)
        MARK_DUPLICATES           (SAMTOOLS_SORT.out.bam, outdir, mark_duplicates_args)
        SAMTOOLS_INDEX            (MARK_DUPLICATES.out.bam, outdir, samtools_index_args)
        BEDTOOLS_GENOMECOV        (MARK_DUPLICATES.out.bam, outdir, bedtools_genomecov_args)
        BEDTOOLS_GENOMECOV_NORM   (MARK_DUPLICATES.out.bam, outdir, bedtools_genomecov_args)



        /* ========================================================================================
            Reports
        ======================================================================================== */

        // Merging channels for MultiQC
        if (!params.skip_qc){

            multiqc_ch = FASTQC.out.report.mix(
                            TRIM_GALORE.out.report.ifEmpty([]),
                            FASTQC2.out.report.ifEmpty([])
                            ).collect()

            if (!params.skip_fastq_screen){
                multiqc_ch = multiqc_ch.mix(
                            FASTQ_SCREEN.out.report.ifEmpty([])
                            ).collect()
            }

            multiqc_ch = multiqc_ch.mix(
                            BOWTIE2.out.stats.ifEmpty([])
                            ).collect()

        } else {

            multiqc_ch = BOWTIE2.out.stats.ifEmpty([]).collect()

        }

        multiqc_ch = multiqc_ch.mix(
            MARK_DUPLICATES.out.metrics.ifEmpty([])
            ).collect()


        MULTIQC (multiqc_ch, outdir, multiqc_args)
}
