#!/bin/bash
# eggd_generate_bed 1.0.0

set -euxo pipefail

main() {

    # download all input files and move to home
    dx-download-all-inputs
    find ~/in -type f -name "*" -print0 | xargs -0 -I {} mv {} ./

    # install python packages from included wheels
    pip install packages/pytz-* packages/numpy-* packages/pandas-*

    # check either panel name or manifest and sample file given
    if [ -z ${panel+x} ]; then
        if [ -z ${manifest_name+x} ] || [ -z ${sample_file_name+x} ]; then
            echo "No panel or manifest and sample identifying file found. Exiting now."
            exit 1;
        else
        # manifest and sample file given, get panel from manifest by sample name
        # get the same ID from the first part of the file name
        sample=$(grep -oP "^[a-zA-Z0-9]*" <<< "$sample_file_name")
        # get the panel(s) from the sample entry in the manifest
        panel=$(grep -w $sample ~/"$manifest_name" | cut -f 2 | sort | uniq | awk '{print}' ORS=', ')

        echo "Sample ID used: $sample"
        fi
    fi

    echo "Using panel(s): $panel"

    # run in empty out dir incase a .bed file is used as input for uploading output bed file
    mkdir ./out && cd ./out

    # generate bed file    
    python3 ~/generate_bed.py -p "$panel" -e ~/"$exons_nirvana_name" -g ~/"$gene_panels_name" -t ~/"$nirvana_genes2transcripts_name"

    bed_file=$(find . -name "*37*.bed" -o -name "*38*.bed")

    echo "Done, uploading BED file: $bed_file"

    output_file=$(dx upload "$bed_file" --brief)

    dx-jobutil-add-output bed_file "$output_file" --class=file

}
