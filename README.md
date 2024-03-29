# eggd_generate_bed (DNAnexus Platform App)

## What does this app do?

Generate a BED file for given panel(s) and / or genes(s).

## What are typical use cases for this app?

This app may be executed as a standalone app.
This app was created for requiring bed files for panels to generate coverage reports from.

## What data are required for this app to run?

This app requires the following files:
- genepanels (from 001_Reference)
- genes2transcripts (from 001_Reference)
- exons (from 001_Reference)

AND also either:

- the panel(s) and / or genes(s) to generate the bed for as a semi-colon ";" separated string, with
    genes being prefixed with an underscore (i.e. DDG2P;_HGNC:4834)

OR

- bioinformatics manifest (from 001_Reference)

- a sample identifying file (this is used to get the appropriate panel from the manifest for the sample,
    this should be named X10000_ as the underscore is used to split and identify the sample). Any sample identifying file may be used (preferentially smaller to reduce download time).

## Optional inputs:
- an "additional regions" file may be provided for a gene or panel. A .tsv file with column headers `chromosome, start, end, gene_panel, transcript` may be provided, where the `gene_panel` column should contain a HGNC gene ID or clinical indication and the `transcript` column should contain a description of the corresponding region. The regions specified in this file will be included in the outputted bed file. Additional columns are allowed with unique headers.
- an integer "flank" which is applied to the start and end of regions in bed file
- an "output_file_prefix" which will prefix the name of the output file

## What does this app output?

This app outputs a BED file.

### This app was made by EMEE GLH
