# generate_bed (DNAnexus Platform App)

## What does this app do?

Generate a BED file for given panel(s) and / or genes(s).

## What are typical use cases for this app?

This app may be executed as a standalone app.
This app was created for requiring bed files for panels to generate coverage reports from.

## What data are required for this app to run?

This app requires the following files:
- genepanels (from 001_Reference)
- genes2transcripts (from 001_Reference)
- exons_nirvana (from 001_Reference) 

AND also either:

- the panel(s) and / or genes(s) to generate the bed for as a comma seperated string, with
    genes being prefixed with an underscore (i.e. DDG2P, _HGNC:4834)

OR

- bioinformatics manifest (from 001_Reference)

- a sample identifying file (this is used to get the appropriate panel from the manifest for the sample, 
    this should be named X10000_ as the underscore is used to split and identify the sample). Any sample identifying file may be used (preferentially smaller to reduce download time).

Optionally, additional regions may be provided for a gene or panel. A tsv file with columns `chrom, start, end, custom string, custom int, gene/panel`. The last column should contain a HGNC gene ID or clinical indication which is used to select regions from this file for generating the bed file. Custom string may be a transcript ID and custom int may be an exon number, or equivalent.

## What does this app output?

This app outputs a BED file.

### This app was made by EMEE GLH
