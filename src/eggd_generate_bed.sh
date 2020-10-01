#!/bin/bash
# eggd_generate_bed 1.0.0

main() {

    # download all input files and move to home
    dx-download-all-inputs
    find ~/in -type f -name "*" -print0 | xargs -0 -I {} mv {} ./

    # check either panel name or manifest and sample file given
    if [[ -z $panel ]] || [[ -z $manifest_name ]] && [[ -z $sample_file_name ]]; then
        echo "No panel or manifest and sample identifying file found. Exiting now."
        exit 1
    fi

    # check if panel string given, if not get from manifest and sample name
    if [[ -z $panel ]]; then
        # get the same ID from the first part of the file name
        # sample=$(echo $sample_file_name | cut -d "_" -f 1)
        sample=$(grep -oP "^[a-zA-Z0-9]*" <<< $sample_file_name)
        # get the panel(s) from the sample entry in the manifest
        panel=$(grep -w $sample ~/$manifest_name | cut -f 2 | sort | uniq | awk '{print}' ORS=', ')

        echo "Sample ID used: " $sample
        echo "Using panel(s) from manifest: " $panel
    fi

    pip3 install pandas==0.24.2

    # run in empty out dir incase a .bed file is used as input for uploading output bed file
    mkdir ./out && cd ./out

    # generate bed file    
    python3 ~/generate_bed.py -p "$panel" -e ~/$exons_nirvana_name -g ~/$gemini_panels_name -t ~/$nirvana_genes2transcripts_name

    echo "Done, uploading BED file."

    bed_file=./*bed

    output_file=$(dx upload $bed_file --brief)

    dx-jobutil-add-output bed_file "$output_file" --class=file

}
