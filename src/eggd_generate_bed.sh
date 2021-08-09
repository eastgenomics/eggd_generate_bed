#!/bin/bash
# eggd_generate_bed

set -exo pipefail

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

            # Check if the sample is present in the manifest
            if ! grep -q $sample ~/"$manifest_name"; then
                echo "Sample ${sample} was not found in the manifest ~/${manifest_name}"
                exit 1
            fi

            # get the panel(s) from the sample entry in the manifest
            panel=$(grep -w $sample ~/"$manifest_name" | cut -f 2 | sort | uniq | awk '{print}' ORS='; ')

            echo "Sample ID used: $sample"
        fi
    fi

    echo "Using panel(s): $panel"

    # run in empty out dir incase a .bed file is used as input for uploading output bed file
    mkdir ./out && cd ./out

    # build optional args string for output prefix and flank
    optional_args=""

    if [ ! -z ${output_file_prefix+x} ]; then
        # replace spaces with underscores if present
        output_file_prefix=$(echo $output_file_prefix | sed 's/ /_/g')
        optional_args+=" -o $output_file_prefix"
    fi

    if [ ! -z ${flank+x} ]; then optional_args+=" -f $flank"; fi

    # generate bed file
    if [ ! -z ${optional_args+x} ]; then
        python3 ~/generate_bed.py -p "$panel" -e ~/"$exons_nirvana_name" -g ~/"$gene_panels_name" -t ~/"$nirvana_genes2transcripts_name" $optional_args
    else
        python3 ~/generate_bed.py -p "$panel" -e ~/"$exons_nirvana_name" -g ~/"$gene_panels_name" -t ~/"$nirvana_genes2transcripts_name"
    fi

    bed_file=$(find . -name "*37*.bed" -o -name "*38*.bed")

    # check if bed file is empty, exit if so
    if [ ! -s $bed_file]; then
        echo "empty bed file generated, exiting now."
        dx-jobutil-report-error "Error: empty bed file generated"
        exit 1
    fi

    echo "Done, uploading BED file: $bed_file"

    output_file=$(dx upload "$bed_file" --brief)

    dx-jobutil-add-output bed_file "$output_file" --class=file

}
