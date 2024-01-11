"""Tests for functions in generate_bed.py"""
import os
import sys
import pytest

sys.path.append(os.path.abspath(
    os.path.join(os.path.realpath(__file__), '../../')))

import generate_bed as gb

TEST_DATA_DIR = (
    os.path.join(os.path.dirname(__file__), 'test_data')
)

@pytest.fixture(name="setup_gene_panels")
def read_in_gene_panels():
    """
    Testing utility to mock output of reading in of test gene_panels file
    Returns:
        pd.Dataframe: df of test_genepanels 
    """
    test_gene_panels_file = f"{TEST_DATA_DIR}/test_genepanels.tsv"
    test_gene_panels = gb.read_to_df(
            test_gene_panels_file, "\t", ["clinical_ind", "panel", "gene"],
            case_change={"column": "clinical_ind", "case": "lower"}
    )
    return test_gene_panels

@pytest.fixture(name="setup_g2t")
def read_in_g2t():
    """
    Testing utility to mock output of reading in of test g2t file
    Returns:
        pd.Dataframe: df of test g2t file
    """
    test_g2t_file = f"{TEST_DATA_DIR}/test_g2t.tsv"
    test_g2t = gb.read_to_df(test_g2t_file, "\t", ["gene", "transcript",
    "clinical_tx","canonical"], case_change={
    "column": "gene", "case": "upper"}
    )
    return test_g2t

@pytest.fixture(name="setup_exons")
def read_in_exons():
    """
    Testing utility to mock output of reading in of test_exons file
    Returns:
        pd.Dataframe: df of test g2t file
    """
    test_exons_file = f"{TEST_DATA_DIR}/test_exons.tsv"
    test_exons = gb.read_to_df(
        test_exons_file, "\t", ["chromosome", "start", "end", "gene", "transcript",
        "exon"]
    )
    return test_exons

class TestReadToDf:
    """Methods to test read_to_df() function from generate_bed.py"""
    def test_gene_panels_read_in_and_case_change(self,setup_gene_panels):
        """
        Method to test if the gene_panels file is read in and cased characters
        in the specified column are changed to the specifed case" 
        """
        assert all(
            clind.islower() for clind in setup_gene_panels["clinical_ind"]), (
                "Not all cased characters were changed to specified case"
        )
    def test_g2t_read_in_case_change(self,setup_g2t):
        """
        Method to test if the gene_panels file is read in and cased characters
        in the specified column are changed to the specifed case
        """
        assert all(gene.isupper() for gene in setup_g2t["gene"]), (
            "Not all cased characters where changed to specified case"
        )

class TestGetGBuild:
    """Methods to test get_g_build() function from generate_bed.py"""
    def test_get_37_genome_build(self):
        """
        Method to test that the correct genome build/file suffix is inferred
        from GRCh37 filenames
        """
        files_37 = ["GCF_000001405.25_GRCh37.p13_genomic.exon_5bp_v2.0.0.tsv",
                    "38_GRCh37_file.tsv", "37_GRCh37_file.tsv"]

        assert all([gb.get_g_build(fn) == "_b37.bed" for fn in files_37]), (
            "Genome buiild was not correctly inferred from a filename"
        )

    def test_get_38_genome_build(self):
        """
        Method to test that the correct genome build/file suffix is inferred
        from GRCh38 filenames
        """
        files_38 = ["GCF_000001405.39_GRCh38.p13_genomic.exon_5bp.tsv",
                    "37_GRCh38_file.tsv", "38_GRCh38_file.tsv"]

        assert all([gb.get_g_build(fn) == "_b38.bed" for fn in files_38]), (
            "Genome buiild was not correctly inferred from a filename"
        )

    def test_no_genome_build(self):
        """
        Method to test that a ValueError is raised if an invalid exons file 
        name is passed" 
        """
        for name in ["build_thirty_seven.tsv","b38_file", "GRCH37.tsv",
                     "grch38.tsv", "G_R_C_h_37.tsv"]:
            with pytest.raises(ValueError):
                gb.get_g_build(name)

    def test_ambiguous_genome_build(self):
        """
        Method to test that an ValueError is raised if an ambiguous exons
        file name is passed
        """
        ambig_file = "GCF_000001405.25_GRCh37_GRCh38.p13_genomic.exon_5bp.tsv"
        with pytest.raises(ValueError):
            gb.get_g_build(ambig_file)

class TestGenesAndReadPanels:
    """
    Methods to test the read_panels() function from generate_bed.py
    """
    def test_gene_not_in_g2t(self, setup_g2t, setup_gene_panels):
        """
        Method to test that the correct AssertionError is raised if a gene 
        provided as part of the "panel_list" input is not present in the g2t df.
        """
        test_gene = "_HGNC:ID_NOT_IN_G2T"

        with pytest.raises(
            AssertionError, match= "not present in genes2transcripts file"
        ):
            gb.read_genes_and_panels(
                panel_list=test_gene, g2t=setup_g2t,
                gene_panels=setup_gene_panels
            )

    def test_panel_not_in_gene_panels(
        self, setup_g2t, setup_gene_panels
    ):
        """
        Method to test that the correct AssertionError is raised if a panel name 
        provided as part of the "panel_list" input is not present in the g2t df.
        """
        test_panel = "NOT A PANEL"
        with pytest.raises(
            AssertionError, match= "not present in gene panels file"
        ):
            gb.read_genes_and_panels(
                panel_list=test_panel, g2t=setup_g2t,
                gene_panels=setup_gene_panels
            )

def test_add_regions_required_headers_present():
    """
    Tests that an AssetionError is raised by read_add_regions() when the required
    headers are missing from the additional regions file
    """
    test_add_regions_file = f"{TEST_DATA_DIR}/test_add_regions.tsv"
    with pytest.raises(AssertionError):
        gb.read_add_regions_file(test_add_regions_file)

class TestGetTranscripts:
    """
    Methods to test the get_transcripts() function from generate_bed.py
    """
    def test_all_genes_have_clinical_transcripts(
        self, setup_g2t, setup_exons
    ):
        """
        Method to test that an AssertionError is raised if not all genes have a
        corresponding clinical transcript and that the error message correctly
        reports which genes have failed   
        """
        # Genes "HGNC:14825" and "HGNC:28208" do not have corresponding
        # transcripts in the test_g2t file whereas "HGNC:4053" does
        test_genes = ["HGNC:14825", "HGNC:28208", "HGNC:4053"]

        with pytest.raises(AssertionError) as e_info:
            gb.get_transcripts(
                g2t=setup_g2t, genes=test_genes,
                exons=setup_exons
            )
        assert ("HGNC:4053" not in str(e_info.value) and
            "HGNC:14825" in str(e_info.value) and 
            "HGNC:28208" in str(e_info.value)), "Incorrect error message given"

    def test_transcript_missing_from_exons(
        self, setup_g2t, setup_exons
    ):
        """
        Method to test that an AssertionError is raised if not all clinical
        transcripts are present in the exons file and that the error message
        correctly reports which transcripts have failed   
        """
        # Genes "HGNC:4053" and "HGNC:329" are present in test_g2t and have
        # corresponding clinical transcripts in test_exons
        # (NM_005101.4 and NM_198576.4, respectively). "HGNC:11918" is present
        # in test_g2t and has a crresponding  clinical transcript (NM_003327.4)
        # however this transcript is missing from test_exons.
        test_genes = ["HGNC:11918", "HGNC:4053", "HGNC:329"]

        with pytest.raises(
            AssertionError, match="NM_003327.4 missing from exons file"
        ):
            gb.get_transcripts(
                g2t=setup_g2t, genes=test_genes,
                exons=setup_exons
            )
