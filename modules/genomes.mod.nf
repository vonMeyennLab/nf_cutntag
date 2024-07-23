#!/usr/bin/env nextflow
nextflow.enable.dsl=2


/* ========================================================================================
    PARAMETERS
======================================================================================== */
params.custom_genome_file = ''


/* ========================================================================================
    FUNCTIONS
======================================================================================== */

// getGenome
def getGenome(name) {

    // Find a file with the same name as the genome in our genomes directory
    scriptDir = workflow.projectDir

    if (params.custom_genome_file){
        
        fileName = params.custom_genome_file

    } else {

        fileName = scriptDir.toString() + "/genomes/" + name + ".genome"

    }

    // die gracefully if the user specified an incorrect genome
    def testFile = new File(fileName)
    if (!testFile.exists()) {
        println("\nFile >>$fileName<< does not exist. Listing available genomes...\n")
        listGenomes()
    }
    else {
        // println ("File $fileName exists.")
    }

    genomeFH = new File (fileName).newInputStream()

    genomeValues = [:]  // initialising map. name is also part of each .genome file

    genomeFH.eachLine {
        sections =  it.split("\\s+",2)
        genomeValues[sections[0]] = sections[1]
    }

    return genomeValues

}
