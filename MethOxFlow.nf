nextflow.enable.dsl=2

/*
* MethOXFlow (MOF) is a pipeline to calculate hydroxymethylation (5hmC) and methylation (5mC) frequencies
* Oxidative Bisulfate (OXBS) and Bisulfate (BS) reads are aligned to a reference to determine 5hmC and 5mC
* Pipeline was largely developed from steps in the Smith Lab MethPipe Manual
* Testing dataset includes OXBS (SRR13665287) and BS (SRR13665293) DNA reads
* extracted from Purkinje cell in Mus musculus
*/


/*
* both reference file and fastq file are paramerized
* if not unputted, they will default to the test mouse data
* mouse data will then only be downloaded if it is not already present
* for testing puposes, reference file, index and fastq files will be provided in project directory
* if needed change parameter directories to local location
*/

// fasta file location for reference file (Test dataset is using mouse reference)
params.reference ='/home/jupyter-andrew.galbraith1/MBB659_Pipeline/Data/Reference/mm10.fa'
// abismal generated index file for reference
params.reference_index ='/home/jupyter-andrew.galbraith1/MBB659_Pipeline/Data/Reference/mm10.index'
// path to all OXBS fastq files for sample
params.BS_fastq = '/home/jupyter-andrew.galbraith1/MBB659_Pipeline/Data/Fastq_subset/BS/*fastq*'
// path to all BS fastq files for sample
params.OXBS_fastq = '/home/jupyter-andrew.galbraith1/MBB659_Pipeline/Data/Fastq_subset/OXBS/*fastq*'
// directory to store pipeline outputs (i.e., meth count files and 5mC/5hmC bed file)
params.output_directory = '/home/jupyter-andrew.galbraith1/MBB659_Pipeline/Out'
// random run boolean
// if true the pipeline is run on a random subselection of fastq files
// for testing purposes, keep false and run on fixed subset data instead
params.RandomRun = false
// number of OXBS and BS fastq files that will be randomly selected
params.TestFastqSize = 3

// create reference file and index variable
reference_file = file(params.reference)
reference_index = file(params.reference_index)

// create channel for all OXBS and BS fastqs
Channel
    .fromPath(params.BS_fastq, checkIfExists: true)
    .set {BS_ch}
Channel 
    .fromPath(params.OXBS_fastq, checkIfExists: true)
    .set {OXBS_ch}

// process to download mouse reference genome
// for testing purposes, skip this step and use reference already installed in directory 
process download_ref { 

    // get and extract mouse reference fasta genome
    output:
        path("${projectDir}/Data/Reference/mm10.fa")
    script:
        """
        wget -P ${projectDir}/Data/Reference http://hgdownload.soe.ucsc.edu/goldenPath/mm10/bigZips/mm10.fa.gz 
        gunzip ${projectDir}/Data/Reference/mm10.fa.gz
        """
}


// process to download test OXBS and BS fastq data
// for testing purposes, skip this step and use data subset already installed in directory 
process download_data {

    // get fastq files from Sequence Read Archive from SRR13665293 (BS) and SRR13665287 (OXBS)
    // this takes a long time to run
    // use data already in project directory to avoid this step
    script:
        """
        fastq-dump --clip -O ${projectDir}/Data/Fastq/BS --read-filter pass SRR13665293
        fastq-dump --clip -O ${projectDir}/Data/Fastq/OXBS --read-filter pass SRR13665287
        """

}

// splits downloaded fastq files to enable pipeline testing on a subset of the data
// this step can be skipped by using subset already available in directory
process create_test_dataset {

    // split fastq file into ones with only 100,000 reads as they will be easier to test
    input:
        path BS
        path OXBS
    
    script:
        """
        split -l 400000 --numeric-suffixes $BS ${projectDir}/Data/Fastq/BS/BS_fastq_
        split -l 400000 --numeric-suffixes $OXBS ${projectDir}/Data/Fastq/OXBS/OXBS_fastq_
        """
}

// create abismal index from fasta reference
process index_reference {
    
    input:
        path reference
        path index

    output:
        path(index)
    
    script: 
        """
        abismalidx $reference $index
        """
}

// seperate fastq files into two files for each read pair
// this is needed for downstream analysis
process get_pair_reads {
    input:
        path OXBS
        path BS
    
    output:
        tuple path("${OXBS}_1.fastq"), path("${OXBS}_2.fastq"), emit: OXBS_paired
        tuple path("${BS}_1.fastq"), path("${BS}_2.fastq"), emit: BS_paired
    
    script:
        """
        sed -ne '1~8{N;N;N;p}' $OXBS > ${OXBS}_1.fastq
        sed -ne '5~8{N;N;N;p}' $OXBS > ${OXBS}_2.fastq
        sed -ne '1~8{N;N;N;p}' $BS > ${BS}_1.fastq
        sed -ne '5~8{N;N;N;p}' $BS > ${BS}_2.fastq
        """
}

// remove potential adaptor sequences from reads
process trim_reads {

    input:
        tuple path(OXBS1), path(OXBS2)
        tuple path(BS1), path(BS2)
    
    output:
        tuple path("$OXBS1"), path("$OXBS2"), emit: OXBS_trim
        tuple path("$BS1"), path("$BS2"), emit: BS_trim
    
    script:
        """
        trim_galore --paired -q 0 --length 0 $OXBS1 $OXBS2
        trim_galore --paired -q 0 --length 0 $BS1 $BS2
        """
}

// align paired reads to abismal reference index
process align_paired_reads {
    input:
        tuple path(OXBS1), path(OXBS2)
        tuple path(BS1), path(BS2)
        path(OXBS_path)
        path(BS_path)
        path(index)
    
    output:
        path "${OXBS_path}.sam", emit: OXBS_sam
        path "${BS_path}.sam", emit: BS_sam
 
    
    script:
        """
        abismal -i $index -o ${OXBS_path}.sam $OXBS1 $OXBS2
        abismal -i $index -o ${BS_path}.sam $BS1 $BS2
        """
    
}

// create bam file from sam alignment file
process create_bam_files {
    input:
        path OXBS_sam
        path BS_sam
    
    output:
        path "${OXBS_sam}.bam", emit: OXBS_bam
        path "${BS_sam}.bam", emit: BS_bam
    
    script:
        """
        samtools view -S -b $OXBS_sam > ${OXBS_sam}.bam
        samtools view -S -b $BS_sam > ${BS_sam}.bam
        """
}

// merge all OXBS and BS files into one each
// aside from testing purposes, splitting fastq and then merging bams enables parallization of alignment
process merge_bam_files {
    input:
        path(OXBS_bam)
        path(BS_bam)
    
    output:
        path("OXBS_all.sam"), emit: OXBS_all_sam
        path("BS_all.sam"), emit: BS_all_sam
    
    script:
       
        """
        samtools merge $OXBS_bam -O SAM -o OXBS_all.sam
        samtools merge $BS_bam -O SAM -o BS_all.sam
        """

}

// format sam file for downstream purposes
process format_sam_file {

    input:
        path(OXBS_sam)
        path(BS_sam)
    
    
    output:
        path("f_${OXBS_sam}"), emit: OXBS_format_sam
        path("f_${BS_sam}"), emit: BS_format_sam
    
    script:
        """
        format_reads -o f_${OXBS_sam} -f abismal $OXBS_sam
        format_reads -o f_${BS_sam} -f abismal $BS_sam
        """
    

}


// sort sam file
process sort_sam {

    input:
        path(OXBS_sam)
        path(BS_sam)
    
    
    output:
        path("s${OXBS_sam}"), emit: OXBS_sorted_sam
        path("s${BS_sam}"), emit: BS_sorted_sam
    
    script:
        """
        samtools sort -O sam -o s${OXBS_sam} $OXBS_sam
        samtools sort -O sam -o s${BS_sam} $BS_sam
        """
    

}

// remove any read duplicates
// these are likely due to over PCR amplification 
process remove_duplicates {

    input:
        path(OXBS_sam)
        path(BS_sam)
    
    
    output:
        path("d${OXBS_sam}"), emit: OXBS_no_dup_sam
        path("d${BS_sam}"), emit: BS_no_dup_sam
    
    script:
        """
        duplicate-remover -S ${OXBS_sam}.txt $OXBS_sam d${OXBS_sam}
        duplicate-remover -S ${BS_sam}.txt $BS_sam d${BS_sam}
        """
    

}

// Get methylation counts from OXBS and BS formatted sam files
// Ouputs to unique meth format see details in methpipe manual
process get_methylation_counts {

    //This step takes awhile regardless of file size
    //Has to write out methylation information to each potential CpG

    input:
        path(OXBS_sam)
        path(BS_sam)
        path(output_directory)
        path(reference)
    
    
    output:
        path("${output_directory}/OXBS.meth"), emit: OXBS_meth
        path("${output_directory}/BS.meth"), emit: BS_meth
    
    script:
        """
        methcounts -c $reference -o ${output_directory}/OXBS.meth $OXBS_sam
        methcounts -c $reference -o ${output_directory}/BS.meth $BS_sam
        """
}

process get_5mC_5hmC {
    input:
        path(OXBS_meth)
        path(BS_meth)
        path(output_directory)
    
    shell:
    // This step takes a few seconds if run on command line but takes awhile when run on pipeline
    // Not sure why this is the case 
    // Only taking CpGs from chromosome 1 to reduce time
    """
        awk '\$1 == "chr1"' $OXBS_meth > ${output_directory}/OXBS_filt.meth
        awk '\$1 == "chr1"' $BS_meth > ${output_directory}/BS_filt.meth
    
        mlml -u ${output_directory}/BS_filt.meth -m ${output_directory}/OXBS_filt.meth -o ${output_directory}/5mC_5hmC.bed
    """
}




workflow {

    // only download reference file if it does not already exist
    // will download to parameterized location
    if (!reference_file.exists()){

        reference_path = download_ref()

    }
    else {
        reference_path = Channel.fromPath(params.reference)
    }

    // only download fastqs if they do not already exist
    if (BS_ch.ifEmpty(false) && !OXBS_ch.ifEmpty(false)){

        download_data()

        // create channels for OXBS and BS fastqs
        Channel
            .fromPath(params.BS_fastq, checkIfExists: true)
            .set {BS_ch}
        Channel 
            .fromPath(params.OXBS_fastq, checkIfExists: true)
            .set {OXBS_ch}
    }

    // random run selects a random set of fastqs to run pipeline on
    // ignore for testing purposes
    if (params.RandomRun){

        // splits fastq files each into ones that only have 100,000 lines
        // these small fastq files are more amenable for testing
        create_test_dataset(BS_ch, OXBS_ch)

        // only subset number of the small fastq files are run based on the TestFastqSize set
        Channel.fromPath('/home/jupyter-andrew.galbraith1/MBB659_Pipeline/Data/Fastq/BS/*fastq_*', checkIfExists: true).randomSample(params.TestFastqSize).set{BS_ch}
        Channel.fromPath('/home/jupyter-andrew.galbraith1/MBB659_Pipeline/Data/Fastq/OXBS/*fastq_*', checkIfExists: true).randomSample(params.TestFastqSize).set{OXBS_ch}
    }

    // index reference if it has not already been done so
    if (!reference_index.exists()){

        reference_index_path = index_reference(params.reference, params.reference_index)

    }
    else{
        reference_index_path = Channel.fromPath(params.reference_index)
    }
    
    // these last 10 steps are the bulk of the pipeline
    // the prior steps were for setting up the project directory
    // by using data already available in the directory, the pipeline will start from here
    // the pipeline from this point on the subset data should take 1-2hr
    get_pair_reads(OXBS_ch, BS_ch)
    OXBS_paired_ch = get_pair_reads.out.OXBS_paired
    BS_paired_ch = get_pair_reads.out.BS_paired
    
    trim_reads(OXBS_paired_ch, BS_paired_ch)
    OXBS_trim_ch = trim_reads.out.OXBS_trim
    BS_trim_ch = trim_reads.out.BS_trim
    
    align_paired_reads(OXBS_trim_ch, BS_trim_ch, OXBS_ch, BS_ch, reference_index_path)
    OXBS_sam_ch = align_paired_reads.out.OXBS_sam
    BS_sam_ch = align_paired_reads.out.BS_sam
    
    create_bam_files(OXBS_sam_ch, BS_sam_ch)
    OXBS_bam_ch = create_bam_files.out.OXBS_bam
    BS_bam_ch = create_bam_files.out.BS_bam
    
    OXBS_bam_collect = OXBS_bam_ch.collect{ " ${it} " }
    BS_bam_collect = BS_bam_ch.collect{ " ${it} " }
    
    merge_bam_files(OXBS_bam_collect, BS_bam_collect)
    OXBS_all_sam_ch = merge_bam_files.out.OXBS_all_sam
    BS_all_sam_ch = merge_bam_files.out.BS_all_sam
    
    format_sam_file(OXBS_all_sam_ch, BS_all_sam_ch)
    OXBS_format_sam_ch = format_sam_file.out.OXBS_format_sam
    BS_format_sam_ch = format_sam_file.out.BS_format_sam
    
    sort_sam(OXBS_format_sam_ch, BS_format_sam_ch)
    OXBS_sort_sam_ch = sort_sam.out.OXBS_sorted_sam
    BS_sort_sam_ch = sort_sam.out.BS_sorted_sam
    
    remove_duplicates(OXBS_sort_sam_ch, BS_sort_sam_ch)
    OXBS_no_dup_sam_ch = remove_duplicates.out.OXBS_no_dup_sam
    BS_no_dup_sam_ch = remove_duplicates.out.BS_no_dup_sam
    
    // here is an output of the paths for both sam files to show the pipeline has worked up to this step
    // after this step the pipeline takes a very long time to run
    OXBS_no_dup_sam_ch.view()
    BS_no_dup_sam_ch.view()
    
    
    get_methylation_counts(OXBS_no_dup_sam_ch, BS_no_dup_sam_ch, params.output_directory, params.reference)
    OXBS_meth_ch = get_methylation_counts.out.OXBS_meth
    BS_meth_ch = get_methylation_counts.out.BS_meth
    
    // this file will contain genomic location follwed by 5mC, 5hmC and unmethylated frequencies for each CpG
    // on the test data the file will be mostly to completely empty
    // this is because there will be little to no overlapping data from OXBS and BS due to the subsampling
    // there is a small error rate in this setp that is not due to the pipeline but rather the methpipe functions
    get_5mC_5hmC(OXBS_meth_ch, BS_meth_ch, params.output_directory)
}