import pandas as pd
import pytest
from pathlib import Path
import os

'''
New test cases requeired, as the existing test cases gives error with new  test function for generate_bed() in test_generate_bed.py script.
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
TEST_DATA_DIR = (
    os.path.join(os.path.dirname(__file__), 'test_data')
)

# create a test exon_df- dtypes for required columns are defined in the generated_bed.py script
@pytest.fixture(name="setup_new_exons")
def make_new_exons():
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
    return f"{TEST_DATA_DIR}/test_exons_generate_bed_v1.3.1.tsv"


# create a test g2t_df- dtypes for required columns are defined in the generated_bed.py script 
@pytest.fixture(name="setup_new_g2t")
def make_new_g2t():
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
            "canonical":str})

    # save dataframe to resources/home/dnanexus/tests/test_data
    df_g2t.to_csv('resources/home/dnanexus/tests/test_data/test_g2t_generate_bed_v1.3.1.tsv', sep="\t",header=False,index=False)
    return f"{TEST_DATA_DIR}/test_g2t_generate_bed_v1.3.1.tsv"


# create a test gene panels file- dtypes for required columns are defined in the generated_bed.py script 
@pytest.fixture(name="setup_new_gene_panels")
def make_new_gene_panels():
    column=["clinical_ind", "panel", "gene"]
    df_gene_panels=pd.DataFrame([
        ["R50.1_Early onset dementia","Early onset dementia","HGNC:10000"],
        ["R50.1_Early onset dementia","Early onset dementia","HGNC:12600"],
        ["R50.1_Early onset dementia","Early onset dementia","HGNC:44506"],
        ["R50.1_Early onset dementia","Early onset dementia","HGNC:93010"]
        ], columns=column).astype(
            {"clinical_ind":str, 
            "panel":str, 
            "gene":str})
    df_gene_panels.to_csv('resources/home/dnanexus/tests/test_data/test_gene_panels_generate_bed_v1.3.1.tsv', sep="\t",header=False,index=False)
    return f"{TEST_DATA_DIR}/test_gene_panels_generate_bed_v1.3.1.tsv"

# create a test additonal regions file- dtypes for required columns are defined in the generated_bed.py script 
@pytest.fixture(name="setup_new_additional_regions")
def make_new_addtional_regions():
    column=["chromosome", "start", "end", "gene_panel","transcript","exons"]
    df_additional_regions=pd.DataFrame([
        ["1",110011,110090,"HGNC:10000","NM_080050.4","1"],
        ["1",100010,100200,"HGNC:12600","NM_010000.2","11"],
        ["2",100001,100190,"HGNC:22506","NM_010006.3","."],
        ["9",289070,289078,"HGNC:96010","NM_181711.1","1"]
        ], columns=column).astype(
            {"chromosome":int, 
                "start":int, 
                "end":int,
                "gene_panel": str,
                "transcript": str,
                "exons": str})
    df_additional_regions.to_csv('resources/home/dnanexus/tests/test_data/test_additional_regions_generate_bed_v1.3.1.tsv', sep="\t",index=False)
    return f"{TEST_DATA_DIR}/test_additional_regions_generate_bed_v1.3.1.tsv"

'''
Make the path locations for expected bed files
'''
@pytest.fixture
def make_expected_bed_path(tmp_path:Path)-> Path:
    # specify nested directory structure 
    expected_bed_path=Path(f"{TEST_DATA_DIR}/expected_beds/")
    print(expected_bed_path)
    #created nested directories
    expected_bed_path.mkdir(parents=True,exist_ok=True)
    return expected_bed_path

'''
This scripts makes expected bed files (with prefix test_expected_bed_file) for 4 tests scenarios 
1. without flank and additional regions
2. with flank
3. with additional regions
4. with flank and additional regions

'''
@pytest.fixture(name="setup_expected_bed")
def make_expected_bed_file(make_expected_bed_path):        
    columns=["chromosome","start","end","transcript"]
    expected_bed_file= pd.DataFrame(
        [
        ["1",110001,110005,"NM_080050.4"],
        ["1",100001,100020,"NM_010000.2"],
        ["1",100001,100190,"NM_010000.3"],
        ["1",289070,289078,"NM_181758.1"]
        ] ,columns=columns
        ).astype({"chromosome": str, "start": int ,"end": int, "transcript": str})

    expected_bed_file.to_csv(f"{make_expected_bed_path}/expected_beds/expected_bed_file.bed", sep="\t", header=False, index=False)
    return f"{make_expected_bed_path}/expected_bed_file.bed"


#make expected_bed_file with flank
@pytest.fixture(name="setup_expected_bed_with_flank")
def make_expected_bed_file_with_flank(make_expected_bed_path):
    # make expected bed file and save as test_expected_bed_file.bed      
    columns=["chromosome","start","end","transcript"]
    expected_bed_file= pd.DataFrame(
        [
        ["1",109601,110405,"NM_080050.4"],
        ["1",99601,100420,"NM_010000.2"],
        ["1",99601,100590,"NM_010000.3"],
        ["1",288670,289478,"NM_181758.1"]
        ] ,columns=columns
        ).astype({"chromosome": str, "start": int ,"end": int, "transcript": str})

    expected_bed_file.to_csv(f"{make_expected_bed_path}/expected_beds/expected_bed_file_with_flank.bed", sep="\t", header=False, index=False)
    return f"{make_expected_bed_path}/expected_bed_file_with_flank.bed"

#make expected_bed_file with additional regions
@pytest.fixture(name="setup_expected_bed_with_additional_regions")
def make_expected_bed_file_with_additional_regions(make_expected_bed_path):
    # make expected bed file and save as test_expected_bed_file.bed    
    columns=["chromosome","start","end","transcript"]
    expected_bed_file= pd.DataFrame(
        [
        ["1",110001,110005,"NM_080050.4"],
        ["1",100001,100020,"NM_010000.2"],
        ["1",100001,100190,"NM_010000.3"],
        ["1",289070,289078,"NM_181758.1"],
        ["1",110011,110090,"NM_080050.4"],
        ["1",100010,100200,"NM_010000.2"]
        ] ,columns=columns
    ).astype({"chromosome": str, "start": int ,"end": int, "transcript": str})
    
    expected_bed_file.to_csv(f"{make_expected_bed_path}/expected_beds/expected_bed_file_with_additional_regions.bed", sep="\t", header=False, index=False)
    return f"{make_expected_bed_path}/expected_bed_file_with_additional_regions.bed"


#make expected_bed_file with additional regions and flank
@pytest.fixture(name="setup_expected_bed_with_additional_regions_flank")
def make_expected_bed_file_with_additional_regions_and_flank(make_expected_bed_path):
    # make expected bed file and save as test_expected_bed_file.bed    
    columns=["chromosome","start","end","transcript"]
    expected_bed_file= pd.DataFrame(
        [
        ["1",109601,110405,"NM_080050.4"],
        ["1",99601,100420,"NM_010000.2"],
        ["1",99601,100590,"NM_010000.3"],
        ["1",288670,289478,"NM_181758.1"],
        ["1",109611,110490,"NM_080050.4"],
        ["1",99610,100600,"NM_010000.2"]
        ] ,columns=columns
    ).astype({"chromosome": str, "start": int ,"end": int, "transcript": str})
    
    expected_bed_file.to_csv(f"{make_expected_bed_path}/expected_beds/expected_bed_file_with_additional_regions_and_flank.bed", sep="\t", header=False, index=False)
    return f"{make_expected_bed_path}/expected_bed_file_with_additional_regions_and_flank.bed"