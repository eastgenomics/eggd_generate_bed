"""
Generates bed file of given panel(s) and/or gene(s) from genepanels,
genes2transcripts and exons file.
"""
import argparse
import pandas as pd


def parse_args():
    """
    Parse arguments given at cmd line.

    Args: None

    Returns:
        - args (Namespace): object containing parsed arguments.
    """

    parser = argparse.ArgumentParser(
        description='Generate BED file from genepanels and exons file.'
    )

    parser.add_argument(
        '-p', '--panel',
        help="Name of panel(s) to generate bed for, for multiple panels these\
            should be provided as a comma separated list. Single genes should\
            be given the same but prefixed with '_'\
            (i.e. --panel panelA, panelB, _gene1, _gene2 etc.)",
        required=True
    )

    parser.add_argument(
        '-e', '--exons', help='exons file', required=True
    )

    parser.add_argument(
        '-g', '--gene_panels', help='gene panels file', required=True
    )

    parser.add_argument(
        '-t', '--g2t', help='genes2transcripts file', required=True
    )

    parser.add_argument(
        '-a', '--additional_regions', default=None,
        help='additional regions file'
    )

    parser.add_argument(
        '-o', '--output', default=None, help='Output file prefix'
    )

    parser.add_argument(
        '-f', '--flank', type=int, default=None,
        help='bp flank to add to each bed file region (optional)'
    )

    args = parser.parse_args()

    return args


def read_to_df(
    file_name: str, sep: str, col_names: list = None, case_change: dict = None,
    required_headers: list = None
) -> pd.DataFrame:
    """
    Read input file (either .tsv or .csv) in as a dataframe and modify a
    column to either all upper or lower case letters.

    Args:
        file_name (str): name of the .csv or .tsv file to be read
        sep (str): separator or delimiter to be used to parse file
        col_names (list): list of column names to be used as headers (optional)
        case_change (dict): dictionary specifying the column and desired
            case to which the column will be changed to upper/lower (optional)

    Raises:
        AssertionError if not all required file headers are present in the
            additional regions file

    Returns:
        pd.DataFrame: df of file.
    """
    dtypes = {
        "name": str, "id": str, "gene": str, "transcript": str,
        "clinical_tx": str, "canonical": str, "chromosome": str,
        "start": int, "end": int
    }

    df = pd.read_csv(file_name, sep=sep, names=col_names, dtype=dtypes)

    if case_change:
        if case_change["case"] == "upper":
            df[case_change["column"]] = df[case_change["column"]].str.upper()
        else:
            df[case_change["column"]] = df[case_change["column"]].str.lower()

    if required_headers:
        assert all([header in df.columns for header in required_headers]), (
            "File doesn't have all the required headers"
        )
    return df


def get_genome_build(exons_file: str) -> str:
    """
    Infer genome build from the exon file's name and return the appropriate
    file suffix to be used.

    Args:
        exons_file (str): exons file filename

    Raises:
        ValueError if ambiguous genome build found in the exon file's name
        ValueError if no genome build is found in the exon file's name

    Returns:
        str: genome build file suffix (either "_37.bed" or "_38.bed")
    """

    if "GRCh37" in exons_file and "GRCh38" in exons_file:
        raise ValueError(f"Ambiguous genome build in exon's file {exons_file}")

    if "GRCh38" in exons_file:
        genome_build = "_b38.bed"
    elif "GRCh37" in exons_file:
        genome_build = "_b37.bed"
    else:
        raise ValueError(
            f"Genome build could not be inferred from {exons_file}"
        )
    return genome_build


def read_genes_and_panels(
    panel_list: str, g2t: pd.DataFrame, gene_panels: pd.DataFrame
):
    """
    Reads in panels/genes and checks if they are in the
    gene_panels/g2t file, respectively.

    Args:
        panel_list (str): semi-colon separated list of panels/genes
        g2t (df): df of genes2transcripts file
        gene_panels (df): df of gene_panels file

    Raises:
        AssertionError if a gene is not present in the genes2transcripts file
        AssertionError if a panel is not present in the gene panels file

    Returns:
        panels (list): list of panels/genes to generate bed for
        genes (list): list of unique genes across all specified panels/genes
    """
    panels = list(filter(None, [x.strip() for x in panel_list.split(";")]))
    genes = []
    for panel in panels:
        if panel.startswith("_"):
            gene = panel.upper().replace("_", "")
            assert gene in g2t["gene"].to_list(), (
                f"Gene {panel} not present in genes2transcripts file"
            )
            genes.append(gene)
        else:
            assert panel.lower() in gene_panels["clinical_ind"].to_list(), (
                f"Panel {panel} not present in gene panels file"
            )
            genes.extend(gene_panels.loc[
                    gene_panels["clinical_ind"] == panel.lower(),
                    "gene"].to_list()
            )

    genes = [x.upper() for x in genes]
    genes = list(set(genes))

    return panels, genes


def get_transcripts(
    g2t: pd.DataFrame, genes: list, exons: pd.DataFrame
) -> list:
    """
    Get transcripts for each gene using g2t.

    Args:
        g2t (pd.DataFrame): df of genes2transcripts file
        genes (list): unique list of genes from all specified panel(s)
        exons (pd.DataFrame): df of exons file

    Raises:
        AssertionError if a selected gene does not have a corresponding
            clinical transcript
        AssertionError if selected transcript is missing from the exons file

    Returns:
        list: list of clinical transcripts corresponding to each gene
    """
    filtered_g2t = g2t.loc[
        (g2t["gene"].isin(genes)) &
        (g2t["clinical_tx"] == "clinical_transcript"), ]

    genes_with_clin_transcripts = set(filtered_g2t["gene"].tolist())

    assert genes_with_clin_transcripts == set(genes), (
        "The following genes do not have corresponding clinical transcripts "
        f"in g2t: {sorted(set(genes).difference(genes_with_clin_transcripts))}"
    )

    # select transcript for each gene in panel genes from entry in g2t
    # get unique in case of duplicates
    transcripts = filtered_g2t["transcript"].unique().tolist()

    # check all selected transcripts in exons file
    for transcript in transcripts:
        assert transcript in exons["transcript"].to_list(), (
            f"{transcript} missing from exons file. Exiting now."
        )

    return transcripts


def generate_bed(
    exons, transcripts, panels, genes, genome_build,
    output_prefix=None, additional_regions=None, flank=None
):
    """
    Generate bed file from transcripts and exons.

    Args:
        - panels (list): name of panel(s) or gene(s)
        - genes (list): unique list of genes from all specified panel(s)
        - transcripts (list): list of transcripts for given set of genes
        - exons (df): df of exons file
        - genome_build (str): file suffix either "_b37.bed" or "_b38.bed"
        - output_prefix (str): Prefix to be added if passed (optional)
        - additional_regions (df) : df of additional_regions file (optional)
        - flank (int) : bp flank to add to each bed file region (optional)

    Returns: None

    Outputs: panel bed file
    """
    # get exons from exons file for transcripts
    # get required columns for bed file
    panel_bed = exons.loc[
        exons["transcript"].isin(transcripts), [
            "chromosome", "start", "end", "transcript"
        ]
    ]

    if additional_regions:
        extra_regions = additional_regions.loc[
            (additional_regions["gene_panel"].isin(panels)) | (
                additional_regions["gene_panel"].isin(genes)),
            ["chromosome", "start", "end", "transcript"]]
        panel_bed = pd.concat([panel_bed, extra_regions], ignore_index=True)

    # apply flank to start and end if given
    if flank:
        print(f"Applying flank of {flank} bp")
        # prevent start becoming -ve where flank is large than distance to 0
        panel_bed.start = panel_bed.start.apply(
            lambda x: x - flank if x - flank >= 0 else 0
        )
        panel_bed.end = panel_bed.end.apply(lambda x: x + flank)

    # write output bed file
    panels = [x.strip(" ").replace(" ", "_") for x in panels]
    panels = [x.replace("/", "-") for x in panels]

    if not output_prefix:
        # Add if statement to minimise long output_prefix as
        # dnanexus cannot handle long filenames
        length_panels = sum(len(i) for i in panels)
        length_dividers = (len(panels) - 1) * 2
        length_output = length_panels + length_dividers
        print("Length of output prefix: " + str(length_output))

        if length_output > 240:
            output_prefix = "".join(panels[0:3])
            output_prefix = output_prefix + "_+" + \
                str(len(panels)-3) + "others"
        else:
            output_prefix = "&&".join(panels)

    if flank:
        # add flank used to output name
        output_prefix = f"{output_prefix}_{flank}bp"

    outfile = output_prefix + genome_build

    panel_bed.to_csv(outfile, sep="\t", header=False, index=False)


def main():
    """
    Main function to generate bed file.
    """
    args = parse_args()
    gene_panels = read_to_df(
        args.gene_panels, "\t", ["clinical_ind", "panel", "gene"],
        case_change={"column": "clinical_ind", "case": "lower"}
    )
    g2t = read_to_df(
        args.g2t, "\t", ["gene", "transcript", "clinical_tx", "canonical"],
        case_change={"column": "gene", "case": "upper"}
    )
    exons = read_to_df(
        args.exons, "\t", ["chromosome", "start", "end", "gene", "transcript",
                           "exon"]
    )
    genome_build = get_genome_build(args.exons)
    panels, genes = read_genes_and_panels(args.panel, g2t, gene_panels)
    transcripts = get_transcripts(g2t, genes, exons)

    if args.additional_regions:
        args.additional_regions = read_to_df(
            file_name=args.additional_regions,
            sep="\t",
            required_headers=["chromosome", "start", "end", "gene_panel",
                              "transcript"]
        )

    generate_bed(
        exons=exons,
        transcripts=transcripts,
        panels=panels,
        genes=genes,
        genome_build=genome_build,
        output_prefix=args.output,
        additional_regions=args.additional_regions,
        flank=args.flank)


if __name__ == "__main__":
    main()
