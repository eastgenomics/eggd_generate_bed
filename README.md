# generate_bed (DNAnexus Platform App)

## What does this app do?

Generate a BED file for given panel(s).

## What are typical use cases for this app?

This app may be executed as a standalone app.
This app was created for requiring bed files for Gemini panels to generate coverage reports from.

## What data are required for this app to run?

This app requires the following files:
- genepanels (from 001_References)
- genes2transcripts (from 001_References)
- exons_nirvana (from 001_References) 

AND also either:

- the panel(s) and / or genes(s) to generate the bed for as a comma seperated string, with
    genes being prefixed with an underscore (i.e. DDG2P, _CFTR)

OR

- bioinformatics manifest (from 001_References)<br>
AND
- a sample identifying file (this is used to get the appropriate panel from the manifest for the sample, 
    this should be named X10000_ as the underscore is used to split and identify the sample). Any sample identifying file may be used (preferentially smaller to reduce download time)



## What does this app output?

This app outputs a BED file.

### This app was made by EMEE GLH
