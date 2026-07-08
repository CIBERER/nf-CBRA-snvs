#!/bin/bash

set -euo pipefail

input=$1        # input.tab.gz
output=$2       # output.tab.gz
maf=$3          # ex. 0.1

zcat "$input" | \
gawk -F'\t' -v OFS='\t' -v maf="$maf" '

BEGIN{
    header_found=0
}

# Mantener las líneas de metadatos
/^##/{
    print
    next
}

# Procesar la cabecera
/^#Uploaded_variation/{

    header_found=1

    for(i=1;i<=NF;i++){
        name=$i
        sub(/^#/,"",name)
        col[name]=i
    }

    use_gnomad = ("gnomADe_AF_grpmax" in col && "gnomADg_AF_grpmax" in col)

    if(use_gnomad){
        eAF   = col["gnomADe_AF_grpmax"]
        eFilt = col["gnomADe_filt"]
        gAF   = col["gnomADg_AF_grpmax"]
        gFilt = col["gnomADg_filt"]
    }else{
        maxAF = col["MAX_AF"]
    }

    print
    next
}

# Ignorar cualquier línea antes de encontrar la cabecera
header_found==0{
    next
}

# Procesar variantes
{

    if(use_gnomad){

        e=$eAF
        g=$gAF

        sub(/,.*/, "", e)
        sub(/,.*/, "", g)

        keepE = (e=="" || e=="." || e=="NA" || e<maf || $eFilt!="PASS")
        keepG = (g=="" || g=="." || g=="NA" || g<maf || $gFilt!="PASS")

        if(keepE && keepG)
            print

    }else{

        x=$maxAF
        sub(/,.*/, "", x)

        if(x=="" || x=="." || x=="NA" || x<maf)
            print
    }
}

' | bgzip -c > "$output"