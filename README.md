# CUT&Tag Sequencing Pipeline

<img width="30%" src="https://raw.githubusercontent.com/nextflow-io/trademark/master/nextflow-logo-bg-light.png" />
<img width="30%" src="https://tower.nf/assets/nf-tower-black.svg" />

A Nextflow pipeline to perform quality control, alignment, and feature coverage of CUT&Tag sequencing data.

>The pipeline was created to run on the [ETH Euler cluster](https://scicomp.ethz.ch/wiki/Euler) and it relies on the server's genome files. Thus, the pipeline needs to be adapted before running it in a different HPC cluster.


## Pipeline steps
1. [FastQC](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/)
2. [FastQ Screen](https://www.bioinformatics.babraham.ac.uk/projects/fastq_screen/)
3. [Trim Galore](https://www.bioinformatics.babraham.ac.uk/projects/trim_galore/)
4. [FastQC](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/)
5. [Bowtie2](https://bowtie-bio.sourceforge.net/bowtie2/)
6. [Samtools sort](https://www.htslib.org/doc/samtools-sort.html)
7. [picard MarkDuplicates](https://broadinstitute.github.io/picard/command-line-overview.html#MarkDuplicates)
8. [Samtools index](https://www.htslib.org/doc/samtools-index.html)
9. [bedtools genomecov](https://bedtools.readthedocs.io/en/latest/content/tools/genomecov.html)
10. [MultiQC](https://multiqc.info/)


## Required parameters

Path to the folder where the FASTQ files are located.
`--input`
``` bash
--input /cluster/work/nme/data/josousa/project/fastq/*fastq.gz
```

Output directory where the files will be saved.
`--outdir`
``` bash
--outdir /cluster/work/nme/data/josousa/project
```


## Input optional parameters

- Option to force the pipeline to assign input as single-end.
`--single_end`

    >_By default, the pipeline detects whether the input files are single-end or paired-end._


## Genomes
- Reference genome used for alignment.
`--genome`

    Available genomes:
    ``` bash
        GRCm39 # Default
        GRCm38
        GRCh38
        GRCh37 
        panTro6
        CHIMP2.1.4
        BDGP6
        susScr11
        Rnor_6.0
        R64-1-1
        TAIR10
        WBcel235
        E_coli_K_12_DH10B
        E_coli_K_12_MG1655
        Vectors
        Lambda
        PhiX
        Mitochondria
    ```

- Option to use a custom genome for alignment by providing an absolute path to a custom genome file.

    ``` bash
    --custom_genome_file '/cluster/work/nme/data/josousa/project/genome/CHM13.genome'
    ```

    Example of a genome file:
    ``` bash    
    name           GRCm39                                                                      
    species        Mouse                                                                       
    fasta          /cluster/work/nme/genomes/Mus_musculus/Ensembl/GRCm39/Sequence/WholeGenomeFasta/           
    bismark        /cluster/work/nme/genomes/Mus_musculus/Ensembl/GRCm39/Sequence/BismarkIndex/               
    bowtie         /cluster/work/nme/genomes/Mus_musculus/Ensembl/GRCm39/Sequence/BowtieIndex/genome          
    bowtie2        /cluster/work/nme/genomes/Mus_musculus/Ensembl/GRCm39/Sequence/Bowtie2Index/genome         
    star           /cluster/work/nme/genomes/Mus_musculus/Ensembl/GRCm39/Sequence/STARIndex/            
    bwa            /cluster/work/nme/genomes/Mus_musculus/Ensembl/GRCm39/Sequence/BWAIndex/genome             
    hisat2         /cluster/work/nme/genomes/Mus_musculus/Ensembl/GRCm39/Sequence/Hisat2Index/genome          
    hisat2_splices /cluster/work/nme/genomes/Mus_musculus/Ensembl/GRCm39/Sequence/Hisat2Index/splice_sites.txt
    gtf            /cluster/work/nme/genomes/Mus_musculus/Ensembl/GRCm39/Annotation/Genes/genes.gtf           
    ```


## FastQ Screen optional parameters

- Option to provide a custom FastQ Screen config file.
    ``` bash
    --fastq_screen_conf '/cluster/work/nme/software/config/fastq_screen.conf' # Default
    ```


## Bowtie2 optional parameters

- Option to suppress SAM records for unaligned reads.
`--bowtie2_no_unal` _Default: true_

- By default, Bowtie2 has the following parameters adapted for CutnTag sequencing:
    ``` bash
    --bowtie2_args="--local --very-sensitive-local --minins 10 --maxins 700"
    ```


## picard MarkDuplicates optional parameters

- Option to not write duplicates to the output file instead of writing them with appropriate flags set.
`--picard_markduplicates_remove_duplicates` _Default: false_

- Option to remove 'optical' duplicates and other duplicates that appear to have arisen from the sequencing process instead of the library preparation process, even if REMOVE_DUPLICATES is false. If REMOVE_DUPLICATES is true, all duplicates are removed and this option is ignored..
`--picard_markduplicates_remove_sequencing_duplicates` _Default: false_


## bedtools genomecov optional parameters

- Option to report depth in BedGraph format, as above (-bg). However with this option, regions with zero coverage are also reported. This allows one to quickly extract all regions of a genome with 0 coverage by applying: 'grep -w 0$' to the output..
`--bedtools_genomecov_bga` _Default: true_


## Skipping options
- Option to skip FastQC, TrimGalore, and FastQ Screen. The first step of the pipeline will be the Bismark alignment. 
`--skip_qc`

- Option to skip FastQ Screen. 
`--skip_fastq_screen`


## Extra arguments
- Option to add extra arguments to [FastQC](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/).
`--fastqc_args`

- Option to add extra arguments to [FastQ Screen](https://www.bioinformatics.babraham.ac.uk/projects/fastq_screen/).
`--fastq_screen_args`

- Option to add extra arguments to [Trim Galore](https://www.bioinformatics.babraham.ac.uk/projects/trim_galore/).
`--trim_galore_args`

- Option to add extra arguments to the [Bowtie2](https://bowtie-bio.sourceforge.net/bowtie2/) aligner.
`--bowtie2_args`

- Option to add extra arguments to [Samtools sort](https://www.htslib.org/doc/samtools-sort.html).
`--samtools_sort_args`

- Option to add extra arguments to [picard MarkDuplicates](https://broadinstitute.github.io/picard/command-line-overview.html#MarkDuplicates).
`--mark_duplicates_args`

- Option to add extra arguments to [Samtools index](https://www.htslib.org/doc/samtools-index.html).
`--samtools_index_args`

- Option to add extra arguments to [bedtools genomecov](https://bedtools.readthedocs.io/en/latest/content/tools/genomecov.html).
`--bedtools_genomecov_args`

- Option to add extra arguments to [MultiQC](https://multiqc.info/).
`--multiqc_args`


## Acknowledgements
This pipeline was adapted from the Nextflow pipelines created by the [Babraham Institute Bioinformatics Group](https://github.com/s-andrews/nextflow_pipelines) and from the [nf-core](https://nf-co.re/) pipelines. We thank all the contributors for both projects. We also thank the [Nextflow community](https://nextflow.slack.com/join) and the [nf-core community](https://nf-co.re/join) for all the help and support.