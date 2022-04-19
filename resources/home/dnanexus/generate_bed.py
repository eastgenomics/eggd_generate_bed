"""
Generates bed file of given panel(s) and/or gene(s) from genepanels,
genes2transcripts and exons_nirvana file.
"""
import argparse
import pandas as pd


def parse_args():
    """
    Pargse arguments given at cmd line.

    Args: None

    Returns:
        - args (Namespace): object containing parsed arguments.
    """

    parser = argparse.ArgumentParser(
        description='Generate BED file from genepanels and exons_nirvana file.'
    )

    parser.add_argument(
        '-p', '--panel',
        help="Name of panel(s) to generate bed for, for multiple panels these\
            should be provided as a comma seperated list. Single genes should\
            be given the same but prefixed with '_'\
            (i.e. --panel panelA, panelB, _gene1, _gene2 etc.)",
        required=True
    )

    parser.add_argument(
        '-e', '--exons_nirvana', help='exons_nirvana file', required=True
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
        '-f', '--flank', type=int,
        help='bp flank to add to each bed file region (optional)'
    )

    args = parser.parse_args()

    return args


def load_files(args):
    """
    Read in files.

    Args:
        - args (Namespace): object containing parsed arguments.

    Returns:
        - panels (list): list of panels to generate bed for
        - gene_panels (df): df of gene_panels file
        - exons_nirvana (df): df of exons_nirvana file
        - g2t (df): df of genes2transcripts file
        - build_38 (bool): check for build of exons nirvana used
    """

    with open(args.gene_panels) as gene_file:
        gene_panels = pd.read_csv(
            gene_file, sep="\t", names=["clinical_ind", "panel", "gene"],
            dtype={"name": str, "id": str, "gene": str},
        )
        gene_panels["clinical_ind"] = gene_panels["clinical_ind"].str.lower()

    with open(args.g2t) as g2t_file:
        g2t = pd.read_csv(
            g2t_file, sep="\t", names=[
                "gene", "transcript", "clinical_tx", "canonical"
            ],
            dtype={
                "gene": str, "transcript": str, "clinical_tx": str,
                "canonical": str
            }
        )
        g2t["gene"] = g2t["gene"].str.upper()

    with open(args.exons_nirvana) as exon_file:
        exons_nirvana = pd.read_csv(
            exon_file, sep="\t",
            names=["chromosome", "start", "end", "gene", "transcript", "exon"],
            dtype={
                "chromosome": str, "start": int, "end": int, "gene": str,
                "transcript": str, "exon": int
            }
        )

    if args.additional_regions:
        with open(args.additional_regions) as extra_regions:
            additional_regions = pd.read_csv(extra_regions, sep="\t")
            headers = ["chromosome", "start", "end", "gene_panel"]
            assert all([x in additional_regions.columns for x in headers]), \
                "Additional regions file doesn't have all the required headers"
    else:
        additional_regions = pd.DataFrame(
            columns=["chromosome", "start", "end", "gene_panel"]
        )

    # check if exons nirvana 37 or 38 used to name output bed
    if "38" in args.exons_nirvana:
        build38 = True
    else:
        build38 = False

    # build list of panels from given string
    panels = args.panel
    panels = list(filter(None, [x.strip() for x in panels.split(";")]))

    # check passed genes in g2t, panels in gene panels
    for panel in panels:
        if panel.startswith("_"):
            assert panel.upper().replace("_", "") in g2t["gene"].to_list(), """
                Gene {} not present in genes2transcripts file""".format(panel)
        else:
            assert panel.lower() in gene_panels["clinical_ind"].to_list(), """\
                Panel {} not present in gene panels file""".format(panel)

    return panels, gene_panels, g2t, exons_nirvana, additional_regions, build38


def generate_bed(
    panels, gene_panels, g2t, exons_nirvana, additional_regions, build38,
    output_prefix, flank=None
):
    """
    Get panel genes from gene_panels for given panel, get transcript
    to use for each gene from g2t then generate bed from transcripts and
    exons_nirvana.

    Args:
        - panels (list): name of panel(s)
        - gene_panels (df): df of gene_panels file
        - exons_nirvana (df): df of exons_nirvana file
        - g2t (df): df of genes2transcripts file
        - build38 (bool): check for build of exons nirvana used
        - output_prefix (str): Prefix to be added if passed

    Returns: None

    Outputs: panel bed file
    """
    # get list of panel genes for each panel
    genes = []

    for panel in panels:
        if panel.startswith("_"):
            # single gene
            genes.append(panel.upper().strip("_"))
        else:
            # panel
            genes.extend(gene_panels.loc[
                    gene_panels["clinical_ind"] == panel.lower(),
                    "gene"].to_list()
            )

    # ensure everything upper case (i.e. instances of lowercase 'orf')
    genes = [x.upper() for x in genes]

    # get unique list of genes across panels
    genes = list(set(genes))

    # select transcript for each gene in panel genes from entry in g2t
    # get unique in case of duplicates
    transcripts = g2t.loc[
            (g2t["gene"].isin(genes)) &
            (g2t["clinical_tx"] == "clinical_transcript"), "transcript"
        ].unique().tolist()

    # check all selected transcripts in exons file
    for transcript in transcripts:
        assert transcript in exons_nirvana["transcript"].to_list(), """\
        {} missing from exons file. Exiting now.""".format(transcript)

    # get exons from exons_nirvana for transcripts
    # get required columns for bed file
    panel_bed = exons_nirvana.loc[
        exons_nirvana["transcript"].isin(transcripts),
            ["chromosome", "start", "end", "transcript"]]

    if len(additional_regions) > 0:
        extra_regions = additional_regions.loc[
            (additional_regions["gene_panel"].isin(panels)) | (
                additional_regions["gene_panel"].isin(genes)),
            ["chromosome", "start", "end", "transcript"]]
        panel_bed = pd.concat([panel_bed, extra_regions], ignore_index=True)


    # apply flank to start and end if given
    if flank:
        print('Applying flank of {} bp'.format(flank))
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
        output_prefix = output_prefix + "_{}bp".format(flank)

    if build38:
        outfile = output_prefix + "_b38.bed"
    else:
        outfile = output_prefix + "_b37.bed"

    panel_bed.to_csv(outfile, sep="\t", header=False, index=False)


def main():
    """
    Main function to generate bed file.
    """
    args = parse_args()

    panels, gene_panels, g2t, exons_nirvana, \
        additional_regions, build38 = load_files(args)

    generate_bed(
        panels, gene_panels, g2t, exons_nirvana,
        additional_regions, build38,
        args.output, args.flank
    )


if __name__ == "__main__":

    main()
