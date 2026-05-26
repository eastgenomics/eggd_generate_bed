'''
This scripts makes new test cases, as the existing test cases gives error with new  test function for generate_bed() in test_generate_bed.py script.
Manually adding new rows to  the existing test case files would also change  formating
For example:
g2t
36  HGNC:10000  NM_080050.4 NOT_CLINICAL_TRANSCRIPT   not_canonical                      NaN            NaN
37      HGNC:12600  NM_010000.2 CLINICAL_TRANSCRIPT       canonical                      NaN            NaN
38  HGNC:44506  NM_010000.3 NOT_CLINICAL_TRANSCRIPT   not_canonical                      NaN            NaN
39      HGNC:93010  NM_181758.1 CLINICAL_TRANSCRIPT       canonical                      NaN            NaN

gene_panels
20  r50.1_early onset dementia  early onset dement...                                                NaN         NaN
21  r50.1_early onset dementia  early onset dement...                                                NaN         NaN
22  r50.1_early onset dementia  early onset dement...                                                NaN         NaN
23  r50.1_early onset dementia  early onset dement...                                                NaN         NaN
'''
import pandas as pd
# create a test exon_df- dtypes for required columns are defined in the generated_bed.py script 
column=["chromosome", "start", "end", "gene", "transcript", "exon"]
exon_df=pd.DataFrame([
         ["1",110001,110005,"HGNC:10000","NM_080050.4",1],
         ["1",100001,100020,"HGNC:12600","NM_010000.2",11],
         ["1",100001,100190,"HGNC:44506","NM_010000.3",6],
         ["1",289070,289078,"HGNC:93010","NM_181758.1",1]],columns=column).astype(
         {"chromosome":str,
                "start":int,
                "end": int,
                "gene": str, 
                "transcript": str,
                "exon" : int})
# save dataframe to resources/home/dnanexus/tests/test_data
exon_df.to_csv('resources/home/dnanexus/tests/test_data/test_exons_generate_bed_v1.3.1.tsv', sep="\t",header=False,index=False)

# create a test g2t_df- dtypes for required columns are defined in the generated_bed.py script 
column=["gene", "transcript", "clinical_tx", "canonical"]
df_g2t=pd.DataFrame([
    ["HGNC:10000","NM_080050.4","clinical_transcript","canonical"],
    ["HGNC:12600","NM_010000.2","clinical_transcript","canonical"],
    ["HGNC:44506","NM_010000.3","clinical_transcript","canonical"],
    ["HGNC:93010","NM_181758.1","clinical_transcript","canonical"]
    ], columns=column).astype(
        {"gene":str, 
         "transcript":str, 
         "clinical_tx":str, 
         "canonical":str}
    )

# save dataframe to resources/home/dnanexus/tests/test_data
df_g2t.to_csv('resources/home/dnanexus/tests/test_data/test_g2t_generate_bed_v1.3.1.tsv', sep="\t",header=False,index=False)

# create a test gene panels file- dtypes for required columns are defined in the generated_bed.py script 
column=["clinical_ind", "panel", "gene"]
df_gene_panels=pd.DataFrame([
        ["R50.1_Early onset dementia","Early onset dementia","HGNC:10000"],
        ["R50.1_Early onset dementia","Early onset dementia","HGNC:12600"],
        ["R50.1_Early onset dementia","Early onset dementia","HGNC:44506"],
        ["R50.1_Early onset dementia","Early onset dementia","HGNC:93010"]
            ]
   , columns=column).astype(
         {"clinical_ind":str, 
          "panel":str, 
          "gene":str}
    )
df_gene_panels.to_csv('resources/home/dnanexus/tests/test_data/test_gene_panels_generate_bed_v1.3.1.tsv', sep="\t",header=False,index=False)

# create a test additonal regions file- dtypes for required columns are defined in the generated_bed.py script 
column=["chromosome", "start", "end", "gene_panel","transcript","exons"]
df_additional_regions=pd.DataFrame([
         ["1",110011,110090,"HGNC:10000","NM_080050.4","1"],
         ["1",100010,100200,"HGNC:12600","NM_010000.2","11"],
         ["2",100001,100190,"HGNC:22506","NM_010006.3","."],
         ["9",289070,289078,"HGNC:96010","NM_181711.1","1"]
            ]
   , columns=column).astype(
         {"chromosome":int, 
          "start":int, 
          "end":int,
          "gene_panel": str,
          "transcript": str,
          "exons": str}
    )
df_additional_regions.to_csv('resources/home/dnanexus/tests/test_data/test_additional_regions_generate_bed_v1.3.1.tsv', sep="\t",index=False)
