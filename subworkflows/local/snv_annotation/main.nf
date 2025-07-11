include { FORMAT2INFO                                                }      from '../../../modules/local/format2info/main'
include { AUTOMAP                                                }      from '../../../modules/local/automap/main'
include { TABIX_TABIX } from '../../../modules/nf-core/tabix/tabix/main'
include { ENSEMBLVEP_VEP                                                }      from '../../../modules/nf-core/ensemblvep/vep/main'
include { VCF_FILTER_MAF } from '../../../modules/local/vcf_filter_maf/main'
include { POSTVEP } from '../../../modules/local/postvep/main'

workflow SNV_ANNOTATION {

    take:
    ch_vcf        // channel (mandatory): [ val(meta), path(vcf), path(tbi) ]
    ch_fasta      // channel (mandatory): [ val(meta), path(fasta) ]
    genome       // channel (mandatory): [ val(genome) ]
    species     // channel (mandatory): [ val(species) ]
    vep_cache_version                      // channel (mandatory): [ val(vep_cache_version) ]
    ch_vep_cache_path                            // channel (mandatory): [ path(cache_path) ]
    ch_vep_custom_extra_files            // channel (optional)  : [ val(meta), path(custom_extra_files) ]
    ch_vep_extra_files                   // channel (optional)  : [ path(extra_files) ]  
    maf
    ch_glowgenes_panel
    ch_glowgenes_sgds

    main:

    ch_versions = Channel.empty()

    ucsc_genome = genome.map { genome_val ->
        switch(genome_val) {
            case 'GRCh38':
                return 'hg38'
            case 'GRCh37':
                return 'hg19'
            default:
                error "Unsupported genome version: ${genome_val}"
        }
    }

    FORMAT2INFO (
         ch_vcf
    )
    ch_versions = ch_versions.mix(FORMAT2INFO.out.versions.first())

    AUTOMAP (
        ch_vcf.map { meta, vcf, tbi -> [meta, vcf] },
        ucsc_genome
    )

    TABIX_TABIX(
        FORMAT2INFO.out.vcf_to_annotate
    )
    ch_versions = ch_versions.mix(TABIX_TABIX.out.versions.first())


    ch_vep_custom_extra_files_and_info = FORMAT2INFO.out.vcf_to_annotate.join(FORMAT2INFO.out.fields).join(TABIX_TABIX.out.tbi).join(ch_vep_custom_extra_files).
        map { tuple ->
        def meta = tuple[0]
        def files = tuple[1..-1].flatten()
        [meta] + files
    }

    //ch_vcf.map { meta, vcf, tbi -> [meta , vcf] }//.view()
    ch_vep = ch_vcf.map { meta, vcf, tbi -> [meta , vcf] }.join(ch_vep_custom_extra_files_and_info).map { items ->
        def meta = items[0]
        def file1 = items[1]
        def restFiles = items[2..-1]  // All files from index 2 to the end
        [
            meta,
            [file1],
            restFiles
        ]
    }
    
    ENSEMBLVEP_VEP (
        ch_vep,
        genome,
        species,
        vep_cache_version,
        ch_vep_cache_path,
        ch_fasta,
        ch_vep_extra_files
    )
    ch_versions = ch_versions.mix(ENSEMBLVEP_VEP.out.versions.first())

    VCF_FILTER_MAF (
        ENSEMBLVEP_VEP.out.tab,
        maf
    )
    
    // Create a complete channel with all metadata from VEP output
    complete_ch = VCF_FILTER_MAF.out.maf_filtered_tab
        .join(AUTOMAP.out.roh_automap_file, remainder: true)
        .map { meta, vep_file, automap_file ->
            // Replace null automap_file with empty list
            automap_file = automap_file ?: []
            [meta, vep_file, automap_file]
        }
    
    //complete_ch.view()
    //ch_glowgenes_panel.view()
    //ch_glowgenes_sgds.view()

    POSTVEP (
        complete_ch,
        maf, 
        ucsc_genome,
        ch_glowgenes_panel,
        ch_glowgenes_sgds
    )

    
    tsv = POSTVEP.out.pvm_tsv

    emit:
    tsv                                         // channel: [ val(meta), path(tsv)]
    versions = ch_versions                     // channel: [ versions.yml ]

}