#!/bin/bash
# eggd_generate_bed 1.0.0

main() {

    echo "Value of panel: '$panel'"
    echo "Value of exons_nirvana: '$exons_nirvana'"
    echo "Value of nirvana_genes2transcripts: '$nirvana_genes2transcripts'"
    echo "Value of gemini_panels: '$gemini_panels'"

    
    dx download "$exons_nirvana" -o exons_nirvana
    dx download "$nirvana_genes2transcripts" -o nirvana_genes2transcripts
    dx download "$gemini_panels" -o gemini_panels

    pip3 install pandas

    python3 generate_bed.py -p $panel -e $exons_nirvana_name -g $gemini_panels_name -t $nirvana_genes2transcripts_name

    bed_file=./*.bed

    bed_file=$(dx upload bed_file --brief)

    dx-jobutil-add-output bed_file "$bed_file" --class=file
}
