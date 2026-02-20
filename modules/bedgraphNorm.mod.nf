#!/usr/bin/env nextflow
nextflow.enable.dsl=2


/* ========================================================================================
	PROCESSES
======================================================================================== */
process NORMALISE_BG {

    label 'normalise'
    tag "$bg" // Adds name to job submission

    input:
        path bg
        val outputdir

    output:
        path "*bedgraph", emit: bedgraph_norm
        publishDir "$outputdir/aligned/bedgraph/normalised/", mode: "link", overwrite: true

    script:
        """
        awk 'FNR==NR{sum+=\$4; next} {printf "%s\t%d\t%d\t%.6f\n",\$1,\$2,\$3, (\$4/sum)*1e6}' ${bg} ${bg} > ${bg}.normalized.bedgraph
                
        for file in *.bedgraph.normalized.bedgraph; do
            mv "\$file" "\${file/.bedgraph.normalized.bedgraph/.normalized.bedgraph}"
        done
        """
}
