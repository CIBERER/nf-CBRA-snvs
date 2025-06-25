include { FORMAT2INFO                                                }      from '../../../modules/local/format2info/main'
include { AUTOMAP                                                }      from '../../../modules/local/automap/main'
include { TABIX_TABIX } from '../../../modules/nf-core/tabix/tabix/main'
include { ENSEMBLVEP_VEP                                                }      from '../../../modules/nf-core/ensemblvep/vep/main'
//include { POSTVEP } from '../../../modules/local/postvep/main'
//include { AUTOMAP } from '../../../modules/local/automap/main' 

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
    genome_ref

    main:

    ch_versions = Channel.empty()
    //ch_vcf.view()
    //ch_extra_files.view()

    FORMAT2INFO (
         ch_vcf
    )
    ch_versions = ch_versions.mix(FORMAT2INFO.out.versions.first())

    AUTOMAP (
        ch_vcf.map { meta, vcf, tbi -> [meta, vcf] },
        genome_ref
    )

    //AUTOMAP.out.roh_automap_file.view()

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
    ch_custom_extra_files_and_info = FORMAT2INFO.out.vcf_to_annotate.join(FORMAT2INFO.out.fields).join(TABIX_TABIX.out.tbi).join(ch_custom_extra_files)//.view()

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
    }.view()
    
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

    // POSTVEP (
    //     ENSEMBLVEP_VEP.out.tab.join(AUTOMAP.out.roh_automap_file),
    //     maf, 
    //     genome
    // )

    
    vcf = ENSEMBLVEP_VEP.out.tab
    vcf.view()

    emit:
    vcf // channel: [ val(meta), path(vcf), path(tbi)]
    versions = ch_versions                     // channel: [ versions.yml ]

}