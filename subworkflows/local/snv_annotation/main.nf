include { FORMAT2INFO                                                }      from '../../../modules/local/format2info/main'
include { AUTOMAP                                                }      from '../../../modules/local/automap/main'
include { TABIX_TABIX } from '../../../modules/nf-core/tabix/tabix/main'
include { ENSEMBLVEP_VEP                                                }      from '../../../modules/nf-core/ensemblvep/vep/main'
include { POSTVEP } from '../../../modules/local/postvep/main'

workflow SNV_ANNOTATION {

    take:
    ch_vcf        // channel (mandatory): [ val(meta), path(vcf), path(tbi) ]
    ch_fasta      // channel (mandatory): [ val(meta), path(fasta) ]
    genome       // channel (mandatory): [ val(genome) ]
    species     // channel (mandatory): [ val(species) ]
    cache_version                      // channel (mandatory): [ val(cache_version) ]
    ch_cache_path                            // channel (mandatory): [ path(cache_path) ]
    ch_custom_extra_files            // channel (optional)  : [ val(meta), path(custom_extra_files) ]
    ch_extra_files                   // channel (optional)  : [ path(extra_files) ]  
    maf

    main:

    ch_versions = Channel.empty()
    //ch_vcf.view()
    //ch_extra_files.view()

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


    // ###################
    // ESTOY INTENTANDO QUE FUNCIONE DESCARGANDO EL CACHE, PERO VA A TARDAR MUCHO, ASÍ QUE PROBAR A PONERLE EL PATH THE VEP.CACHE - '/mnt/genetica5/vep_2024/ 
    // ASSEMBLY hg38
    // genome = 'GRCh38'
    // specie = "homo_sapiens" 

    //cache.view()

    //ch_custom_extra_files.view()
    //ch_custom_extra_files_and_info = ch_custom_extra_files.join(FORMAT2INFO.out.vcf_to_annotate).join(FORMAT2INFO.out.fields).join(TABIX_TABIX.out.tbi)//.view()
    ch_custom_extra_files_and_info = FORMAT2INFO.out.vcf_to_annotate.join(FORMAT2INFO.out.fields).join(TABIX_TABIX.out.tbi).join(ch_custom_extra_files).
        map { tuple ->
        def meta = tuple[0]
        def files = tuple[1..-1].flatten()
        [meta] + files
    }//.view()

    //ch_vcf.map { meta, vcf, tbi -> [meta , vcf] }//.view()
    ch_vep = ch_vcf.map { meta, vcf, tbi -> [meta , vcf] }.join(ch_custom_extra_files_and_info).map { items ->
        def meta = items[0]
        def file1 = items[1]
        def restFiles = items[2..-1]  // All files from index 2 to the end
        [
            meta,
            [file1],
            restFiles
        ]
    }//.view()
    
    ENSEMBLVEP_VEP (
        ch_vep,
        genome,
        species,
        cache_version,
        ch_cache_path,
        ch_fasta,
        ch_extra_files
    )
    ch_versions = ch_versions.mix(ENSEMBLVEP_VEP.out.versions.first())
    
    ENSEMBLVEP_VEP.out.tab.view()

    // Create a complete channel with all metadata from VEP output
    complete_ch = ENSEMBLVEP_VEP.out.tab
        .join(AUTOMAP.out.roh_automap_file, remainder: true)
        .map { meta, vep_file, automap_file ->
            // Replace null automap_file with empty list
            automap_file = automap_file ?: []
            [meta, vep_file, automap_file]
        }
    
    //complete_ch.view()

    POSTVEP (
        complete_ch,
        maf, 
        ucsc_genome
    )

    
    tsv = POSTVEP.out.pvm_tsv
    tsv.view()

    emit:
    tsv // channel: [ val(meta), path(vcf), path(tbi)]
    versions = ch_versions                     // channel: [ versions.yml ]

}