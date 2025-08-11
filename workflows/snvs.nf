/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    PRINT PARAMS SUMMARY
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { paramsSummaryLog; paramsSummaryMap } from 'plugin/nf-schema'

def logo = NfcoreTemplate.logo(workflow, params.monochrome_logs)
def citation = '\n' + WorkflowMain.citation(workflow) + '\n'
def summary_params = paramsSummaryMap(workflow)

// Print parameter summary log to screen
log.info logo + paramsSummaryLog(workflow) + citation

WorkflowSnvs.initialise(params, log)



/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    VALIDATE INPUTS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

// Check mandatory parameters

ch_fasta   = params.fasta ? Channel.fromPath(params.fasta).map{ it -> [ [id:it.baseName], it ] }.collect() : Channel.empty() 
ch_fai     = params.fai ? Channel.fromPath(params.fai).map{ it -> [ [id:it.baseName], it ] }.collect() : Channel.empty()
ch_snps    = params.known_snps ? Channel.fromPath(params.known_snps).collect() : Channel.value([])
ch_snps_tbi = params.known_snps_tbi ? Channel.fromPath(params.known_snps_tbi).collect() : Channel.empty()


//ch_assembly = params.assembly ? Channel.value(params.assembly) : ch_fasta.map { meta, fasta -> meta.id }.first() 
// lo había puesto así porque solo era para poner el id al nombre final del vcf (si no estaba, ponia el nombre del fasta), pero ahora lo he cambiado porque también lo vamos a usar para el VEP
ch_assembly = Channel.value(params.assembly)

//deep variant parameters
ch_gzi = Channel.of([[],[]]).first()
ch_par_bed = params.ch_par_bed ? Channel.fromPath(params.ch_par_bed, checkIfExists: true).map { file -> [ [:], file ] }.collect() : Channel.of([[:], []]).first()


/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    CONFIG FILES
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

ch_multiqc_config          = Channel.fromPath("$projectDir/assets/multiqc_config.yml", checkIfExists: true)
ch_multiqc_custom_config   = params.multiqc_config ? Channel.fromPath( params.multiqc_config, checkIfExists: true ) : Channel.empty()
ch_multiqc_logo            = params.multiqc_logo   ? Channel.fromPath( params.multiqc_logo, checkIfExists: true ) : Channel.empty()
ch_multiqc_custom_methods_description = params.multiqc_methods_description ? file(params.multiqc_methods_description, checkIfExists: true) : file("$projectDir/assets/methods_description_template.yml", checkIfExists: true)

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT LOCAL MODULES/SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

//
// SUBWORKFLOW: Consisting of a mix of local and nf-core/modules
//
include { INPUT_CHECK } from '../subworkflows/local/input_check'
include { MAPPING } from '../subworkflows/local/mapping'
include { GATK_VCF } from '../subworkflows/local/gatk_vcf'
include { DRAGEN_VCF } from '../subworkflows/local/dragen_vcf'
include { VCF_MERGE_VARIANTCALLERS } from '../subworkflows/local/vcf_merge_variantcallers'
include { DEEP_VARIANT_VCF           } from '../subworkflows/local/deep_variant_vcf'
include { SNV_ANNOTATION } from '../subworkflows/local/snv_annotation'
include { GATK_TRIO_VCF } from '../subworkflows/local/gatk_trio_vcf'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT NF-CORE MODULES/SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

//
// MODULE: Installed directly from nf-core/modules
//
include { FASTQC                      } from '../modules/nf-core/fastqc/main'
include { MULTIQC                     } from '../modules/nf-core/multiqc/main'
include { CUSTOM_DUMPSOFTWAREVERSIONS } from '../modules/nf-core/custom/dumpsoftwareversions/main'

include { BWA_INDEX } from '../modules/nf-core/bwa/index/main'
include { PICARD_CREATESEQUENCEDICTIONARY } from '../modules/nf-core/picard/createsequencedictionary/main'
include { GATK4_COMPOSESTRTABLEFILE } from '../modules/nf-core/gatk4/composestrtablefile/main'
include { GATK4_CALIBRATEDRAGSTRMODEL } from '../modules/nf-core/gatk4/calibratedragstrmodel/main'
include { ENSEMBLVEP_DOWNLOAD } from '../modules/nf-core/ensemblvep/download/main'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

// Info required for completion email and summary
def multiqc_report = []

workflow SNVS {

    ch_versions = Channel.empty()

    //
    // SUBWORKFLOW: Read in samplesheet, validate and stage input files
    //
    INPUT_CHECK (
        file(params.input)
    )
    ch_versions = ch_versions.mix(INPUT_CHECK.out.versions)
    // TODO: OPTIONAL, you can use nf-validation plugin to create an input channel from the samplesheet with Channel.fromSamplesheet("input")
    // See the documentation https://nextflow-io.github.io/nf-validation/samplesheets/fromSamplesheet/
    // ! There is currently no tooling to help you write a sample sheet schema
    //INPUT_CHECK.out.reads.view()
    // 
    // MODULE: Run FastQC
    //
    FASTQC (
        INPUT_CHECK.out.reads
    )
    ch_versions = ch_versions.mix(FASTQC.out.versions.first())
     
    if (params.index) { 
        ch_index = Channel.fromPath(params.index).map{ it -> [ [id:it.baseName], it ] }.collect()
    } else { 
        BWA_INDEX (ch_fasta)
        ch_index = BWA_INDEX.out.index
    }

    if (params.refdict) { 
        ch_refdict = Channel.fromPath(params.refdict).map{ it -> [ [id:it.baseName], it ] }.collect()
    } else { 
        PICARD_CREATESEQUENCEDICTIONARY (ch_fasta)
        ch_refdict = PICARD_CREATESEQUENCEDICTIONARY.out.reference_dict
    }

    if (params.reference_str) { 
        ch_ref_str = Channel.fromPath(params.reference_str).collect()
    } else { 
        GATK4_COMPOSESTRTABLEFILE (
            ch_fasta.map {meta, fasta -> [fasta] },
            ch_fai.map {meta, fai -> [fai]  },
            ch_refdict.map {meta, dict -> [dict] }
        )
        ch_ref_str = GATK4_COMPOSESTRTABLEFILE.out.str_table
    }

    // In your main workflow, ensure intervals are created for all samples
    ch_intervals = params.intervals ? 
        INPUT_CHECK.out.reads.map{ meta, fastqs -> tuple(meta, file(params.intervals)) } : 
        INPUT_CHECK.out.reads.map{ meta, fastqs -> tuple(meta, []) }



    MAPPING (
        INPUT_CHECK.out.reads,
        ch_intervals,
        ch_index,
        ch_fasta,
        ch_fai,
        ch_refdict,
        ch_snps,
        ch_snps_tbi
    )
    
    ///////////// TODO: Esto después quitarlo, es solo para probar que funciona el GATK4 de TRIOS ////////////////
    ch_intervals_genomicsdbimport = Channel.fromPath(params.genomicsdbimport_interval).collect()
    //no_intervals = params.intervals ? false : true
    ch_ped = INPUT_CHECK.out.ped.unique()

    if (params.trio_analysis) {
        GATK_TRIO_VCF (
            MAPPING.out.bam,
            ch_intervals,
            ch_fasta,
            ch_fai,
            ch_refdict,
            ch_snps.map{ it -> [ [id:it.baseName], it ] }.collect(),
            ch_snps_tbi.map{ it -> [ [id:it.baseName], it ] }.collect(),
            //Channel.fromList([tuple([ id: 'dbsnp'],[])]).collect(),
            //Channel.fromList([tuple([ id: 'dbsnp_tbi'],[])]).collect(),
            ch_intervals_genomicsdbimport, // ch_intervals_genomicsdbimport
            //no_intervals, // no_intervals
            ch_ped // ch_ped            
        )

        vcf_file = GATK_TRIO_VCF.out.vcf

    } else { ///////////// start OF GATK dragen etc IF BLOCK ////////////////

    GATK_VCF (
        MAPPING.out.bam,
        ch_intervals,
        ch_fasta,
        ch_fai,
        ch_refdict,
        Channel.fromList([tuple([ id: 'dbsnp'],[])]).collect(),
        Channel.fromList([tuple([ id: 'dbsnp_tbi'],[])]).collect()
    )
    
    DEEP_VARIANT_VCF (
        MAPPING.out.bam,
        ch_intervals,
        ch_fasta,
        ch_fai,
        ch_gzi,
        ch_par_bed
    )

   DRAGEN_VCF (
        MAPPING.out.bam, 
        ch_fasta,
        ch_fai,
        ch_refdict,
        GATK4_COMPOSESTRTABLEFILE.out.str_table,
        ch_intervals,
        Channel.fromList([tuple([ id: 'dbsnp'],[])]).collect(),
        Channel.fromList([tuple([ id: 'dbsnp_tbi'],[])]).collect()
    )

    ch_gatk = params.run_gatk ? GATK_VCF.out.vcf : Channel.empty()
    ch_dragstr = params.run_dragen ? DRAGEN_VCF.out.vcf : Channel.empty()
    ch_deepvariant = params.run_deepvariant ? DEEP_VARIANT_VCF.out.vcf : Channel.empty()

    ch_vcfs_for_merge = ch_gatk.join(ch_dragstr).join(ch_deepvariant)

    VCF_MERGE_VARIANTCALLERS (
        ch_vcfs_for_merge,   
        ch_fasta,
        ch_fai,
        ch_intervals,
        ch_assembly
    )

    vcf_file = VCF_MERGE_VARIANTCALLERS.out.vcf

    } //////// END OF GATK_TRIO_VCF IF BLOCK ////////

    ch_custom_extra_files = params.custom_extra_files ? vcf_file.map{ meta, vcf, tbi -> tuple(meta, file(params.custom_extra_files)) } : vcf_file.map{ meta, vcf, tbi -> tuple(meta, []) }
    ch_extra_files = params.extra_files ? Channel.fromPath(params.extra_files, checkIfExists: true).collect() : Channel.value([])

    // Conditionally add files using mix
    if (params.plugins_dir) {
        ch_extra_files = ch_extra_files.mix(Channel.fromPath("${params.plugins_dir}", checkIfExists: true)).collect()
    }

    ch_glowgenes_panel = params.glowgenes_panel ? Channel.fromPath(params.glowgenes_panel, checkIfExists: true).collect() : Channel.value([])
    ch_glowgenes_sgds = params.glowgenes_sgds ? Channel.fromPath(params.glowgenes_sgds, checkIfExists: true).collect() : Channel.value([])

    if (params.vep_cache_path) { ch_vep_cache_path = Channel.fromPath(params.vep_cache_path, checkIfExists: true).collect() } else { 
        // Define your meta_vep
        def meta_vep = [id: "vep_${params.assembly}", assembly: params.assembly]
        if (params.refseq_cache) {
            ch_vep_download = Channel.of([meta_vep, params.assembly, "${params.species}_refseq", params.vep_cache_version])
        } else {
            ch_vep_download = Channel.of([meta_vep, params.assembly, params.species, params.vep_cache_version])
        }
        ENSEMBLVEP_DOWNLOAD (
            ch_vep_download
            )
        ch_vep_cache_path = ENSEMBLVEP_DOWNLOAD.out.cache.map{ meta, cache -> [cache] }.collect()
    }

    ch_vep_cache_version = params.vep_cache_version ? Channel.value(params.vep_cache_version) : Channel.value([])

    SNV_ANNOTATION (
        vcf_file,
        ch_fasta,
        ch_assembly,
        params.species,
        ch_vep_cache_version,
        ch_vep_cache_path,
        ch_custom_extra_files,
        ch_extra_files,
        params.maf,
        ch_glowgenes_panel,
        ch_glowgenes_sgds
    )



    CUSTOM_DUMPSOFTWAREVERSIONS (
        ch_versions.unique().collectFile(name: 'collated_versions.yml')
    )

    //
    // MODULE: MultiQC
    //
    workflow_summary    = WorkflowSnvs.paramsSummaryMultiqc(workflow, summary_params)
    ch_workflow_summary = Channel.value(workflow_summary)

    methods_description    = WorkflowSnvs.methodsDescriptionText(workflow, ch_multiqc_custom_methods_description, params)
    ch_methods_description = Channel.value(methods_description)

    ch_multiqc_files = Channel.empty()
    ch_multiqc_files = ch_multiqc_files.mix(ch_workflow_summary.collectFile(name: 'workflow_summary_mqc.yaml'))
    ch_multiqc_files = ch_multiqc_files.mix(ch_methods_description.collectFile(name: 'methods_description_mqc.yaml'))
    ch_multiqc_files = ch_multiqc_files.mix(CUSTOM_DUMPSOFTWAREVERSIONS.out.mqc_yml.collect())
    ch_multiqc_files = ch_multiqc_files.mix(FASTQC.out.zip.collect{it[1]}.ifEmpty([]))

    MULTIQC (
        ch_multiqc_files.collect(),
        ch_multiqc_config.toList(),
        ch_multiqc_custom_config.toList(),
        ch_multiqc_logo.toList()
    )
    multiqc_report = MULTIQC.out.report.toList()
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    COMPLETION EMAIL AND SUMMARY
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow.onComplete {
    if (params.email || params.email_on_fail) {
        NfcoreTemplate.email(workflow, params, summary_params, projectDir, log, multiqc_report)
    }
    NfcoreTemplate.dump_parameters(workflow, params)
    NfcoreTemplate.summary(workflow, params, log)
    if (params.hook_url) {
        NfcoreTemplate.IM_notification(workflow, params, summary_params, projectDir, log)
    }
}

workflow.onError {
    if (workflow.errorReport.contains("Process requirement exceeds available memory")) {
        println("🛑 Default resources exceed availability 🛑 ")
        println("💡 See here on how to configure pipeline: https://nf-co.re/docs/usage/configuration#tuning-workflow-resources 💡")
    }
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
