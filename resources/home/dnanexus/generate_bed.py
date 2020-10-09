"""
Generates bed file of given panel(s) and/or gene(s) from genepanels,
genes2transcripts and exons_nirvana file.

Jethro Rainford
200928
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
            gene_file, sep="\t", names=["panel", "id", "gene"],
            low_memory=False,
        )
        gene_panels["panel"] = gene_panels["panel"].str.lower()
        gene_panels["gene"] = gene_panels["gene"].str.upper()

    with open(args.exons_nirvana) as exon_file:
        exons_nirvana = pd.read_csv(
            exon_file, sep="\t", low_memory=False,
            names=["chromosome", "start", "end", "gene", "transcript", "exon"]
        )

    # check if exons nirvana 37 or 38 used to name output bed
    if "38" in args.exons_nirvana:
        build38 = True
    else:
        build38 = False

    with open(args.g2t) as g2t_file:
        g2t = pd.read_csv(
            g2t_file, sep="\t", low_memory=False, names=["gene", "transcript"]
        )
        g2t["gene"] = g2t["gene"].str.upper()

    # build list of panels from given string
    panels = args.panel
    panels = list(filter(None, [x.strip() for x in panels.split(",")]))

    # check passed genes in g2t, panels in gene panels
    for panel in panels:
        if "_" in panel:
            assert panel.upper().replace("_", "") in g2t["gene"].to_list(), """
                Gene {} not present in genes2transcripts file""".format(panel)
        else:
            assert panel.lower() in gene_panels["panel"].to_list(), """\
                Panel {} not present in gene panels file""".format(panel)

    return panels, gene_panels, exons_nirvana, g2t, build38


def generate_bed(panels, gene_panels, exons_nirvana, g2t, build38):
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

    Returns: None

    Outputs: panel bed file
    """
    # get list of panel genes for each panel
    genes = []

    for panel in panels:
        if "_" in panel:
            # single gene
            genes.append(panel.upper().strip("_"))
        else:
            genes.extend(
                gene_panels.loc[
                    gene_panels["panel"] == panel.lower()]["gene"].to_list()
            )

    # ensure everything upper case (i.e. instances of lowercase 'orf')
    genes = [x.upper() for x in genes]

    # get unique list of genes across panels
    genes = list(set(genes))

    # get list of transcripts for panel genes
    transcripts = g2t[g2t["gene"].isin(genes)]["transcript"].to_list()

    # get unique in case of duplicates
    transcripts = list(set(transcripts))

    # get exons from exons_nirvana for transcripts
    exons = exons_nirvana[exons_nirvana["transcript"].isin(transcripts)]

    # get required columns for bed file
    panel_bed = exons[["chromosome", "start", "end", "transcript"]]

    # write output bed file
    panels = [x.strip(" ").replace(" ", "_") for x in panels]
    outfile = "&".join(panels)

    if build38:
        outfile = "&&".join(panels) + "_b38.bed"
    else:
        outfile = "&&".join(panels) + "_b37.bed"

    panel_bed.to_csv(outfile, sep="\t", header=False, index=False)


def main():
    """
    Main function to generate bed file.
    """
    args = parse_args()

    panels, gene_panels, exons_nirvana, g2t, build38 = load_files(args)

    generate_bed(panels, gene_panels, exons_nirvana, g2t, build38)


if __name__ == "__main__":

    main()
