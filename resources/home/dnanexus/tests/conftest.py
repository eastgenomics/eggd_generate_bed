import pandas as pd
import pytest
from pathlib import Path
import os

'''
New test cases required, as the existing test cases gives error with new  test function for generate_bed() in test_generate_bed.py script.
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


# create a test exon_df- dtypes for required columns are defined in the generated_bed.py script
@pytest.fixture(name="test_new_exons")
def new_exons(tmp_path: Path) -> str:
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
    # save dataframe to  <class 'pathlib.PosixPath'> which is passed to to_csv (that accepts path-like object)
    out=tmp_path / "test_exons_generate_bed_v1.3.1.tsv"
    exon_df.to_csv(out, sep="\t",header=False,index=False)
    return str(out)


# create a test g2t_df- dtypes for required columns are defined in the generated_bed.py script 
@pytest.fixture(name="test_new_g2t")
def new_g2t(tmp_path: Path) -> str:
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

    # save dataframe  
    out= tmp_path / "test_g2t_generate_bed_v1.3.1.tsv"
    df_g2t.to_csv(out, sep="\t",header=False,index=False)
    return out


# create a test gene panels file- dtypes for required columns are defined in the generated_bed.py script 
@pytest.fixture(name="test_new_gene_panels")
def new_gene_panels(tmp_path: Path) -> str:
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
    out=tmp_path / "test_gene_panels_generate_bed_v1.3.1.tsv"
    df_gene_panels.to_csv(out, sep="\t",header=False,index=False)
    return out

# create a test additonal regions file- dtypes for required columns are defined in the generated_bed.py script 
@pytest.fixture(name="test_new_additional_regions")
def new_additional_regions(tmp_path: Path) -> str:
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
    out=tmp_path / "test_additional_regions_generate_bed_v1.3.1.tsv"
    df_additional_regions.to_csv(out, sep="\t",index=False)
    return out

'''
Make the path locations for expected bed files
'''
@pytest.fixture
def expected_bed_path(tmp_path:Path)-> Path:
    # specify nested directory structure 
    bed_path=tmp_path / "expected_beds"
    print(bed_path)
    #created nested directories
    bed_path.mkdir(parents=True,exist_ok=True)
    return bed_path

'''
This scripts makes expected bed files (with prefix test_expected_bed_file) for 4 tests scenarios 
1. without flank and additional regions
2. with flank
3. with additional regions
4. with flank and additional regions

'''
@pytest.fixture(name="test_expected_bed")
def expected_bed_file(expected_bed_path):        
    columns=["chromosome","start","end","transcript"]
    expected_bed_file= pd.DataFrame(
        [
        ["1",110001,110005,"NM_080050.4"],
        ["1",100001,100020,"NM_010000.2"],
        ["1",100001,100190,"NM_010000.3"],
        ["1",289070,289078,"NM_181758.1"]
        ] ,columns=columns
        ).astype({"chromosome": str, "start": int ,"end": int, "transcript": str})

    expected_bed_file.to_csv(f"{expected_bed_path}/expected_bed_file.bed", sep="\t", header=False, index=False)
    return f"{expected_bed_path}/expected_bed_file.bed"


#make expected_bed_file with flank
@pytest.fixture(name="test_expected_bed_with_flank")
def expected_bed_file_with_flank(expected_bed_path):
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

    expected_bed_file.to_csv(f"{expected_bed_path}/expected_bed_file_with_flank.bed", sep="\t", header=False, index=False)
    return f"{expected_bed_path}/expected_bed_file_with_flank.bed"

#make expected_bed_file with additional regions
@pytest.fixture(name="test_expected_bed_with_additional_regions")
def expected_bed_file_with_additional_regions(expected_bed_path):
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
    
    expected_bed_file.to_csv(f"{expected_bed_path}/expected_bed_file_with_additional_regions.bed", sep="\t", header=False, index=False)
    return f"{expected_bed_path}/expected_bed_file_with_additional_regions.bed"


#make expected_bed_file with additional regions and flank
@pytest.fixture(name="test_expected_bed_with_additional_regions_flank")
def expected_bed_file_with_additional_regions_and_flank(expected_bed_path):
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
    
    expected_bed_file.to_csv(f"{expected_bed_path}/expected_bed_file_with_additional_regions_and_flank.bed", sep="\t", header=False, index=False)
    return f"{expected_bed_path}/expected_bed_file_with_additional_regions_and_flank.bed"