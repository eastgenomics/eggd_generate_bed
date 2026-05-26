import pandas as pd
import os

'''
This scripts makes expected bed files (with prefix test_expected_bed_file) for 4 tests scenarios 
1. without flank and additional regions
2. with flank
3. with additional regions
4. with flank and additional regions

'''
def make_expected_bed_file(test_data_dir):        
    columns=["chromosome","start","end","transcript"]
    expected_bed_file= pd.DataFrame(
        [
        ["1",110001,110005,"NM_080050.4"],
        ["1",100001,100020,"NM_010000.2"],
        ["1",100001,100190,"NM_010000.3"],
        ["1",289070,289078,"NM_181758.1"]
        ] ,columns=columns
        ).astype({"chromosome": str, "start": int ,"end": int, "transcript": str})

    expected_bed_file.to_csv(f"{test_data_dir}/expected_beds/test_expected_bed_file.bed", sep="\t", header=False, index=False)

#make expected_bed_file with flank
def make_expected_bed_file_with_flank(test_data_dir):
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


    expected_bed_file.to_csv(f"{test_data_dir}/expected_beds/test_expected_bed_file_with_flank.bed", sep="\t", header=False, index=False)

#make expected_bed_file with additional regions
def make_expected_bed_file_with_additional_regions(test_data_dir):
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

       
    expected_bed_file.to_csv(f"{test_data_dir}/expected_beds/test_expected_bed_file_with_additional_regions.bed", sep="\t", header=False, index=False)

#make expected_bed_file with additional regions and flank
def make_expected_bed_file_with_additional_regions_and_flank(test_data_dir):
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

       
    expected_bed_file.to_csv(f"{test_data_dir}/expected_beds/test_expected_bed_file_with_additional_regions_and_flank.bed", sep="\t", header=False, index=False)


TEST_DATA_DIR=(
    os.path.join(os.path.dirname(__file__), 'test_data')
)

make_expected_bed_file(test_data_dir=TEST_DATA_DIR)
make_expected_bed_file_with_flank(test_data_dir=TEST_DATA_DIR)
make_expected_bed_file_with_additional_regions(test_data_dir=TEST_DATA_DIR)
make_expected_bed_file_with_additional_regions_and_flank(test_data_dir=TEST_DATA_DIR)
