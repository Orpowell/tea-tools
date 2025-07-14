import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import csv
import logging
import sys


def load_transcript_quantification(quants: str) -> dict[str, float]:
    with open(quants) as file:
        reader = csv.reader(file, delimiter="\t")

        return {row[0]: float(row[1]) for row in reader}


def visualise_transcripts_expression(
    quantification: str,
    graph_out: str,
    genes_of_interest=[str],
    aliases=[str],
    show=False,
    percentage: float = 0.25,
    percentage_colour: str = "ffa500",
    no_legend: bool = False,
):
    if (percentage > 1) or (percentage < 0):
        logging.error(f" percentage ({percentage}) must be between 0 and 1.")
        sys.exit(1)

    if len(genes_of_interest) != len(aliases):
        logging.error(
            f" alias ({len(aliases)}) and transcript-of-interest ({len(genes_of_interest)}) lists are not the same length!"
        )
        logging.info("Shutting down...")
        sys.exit(1)

    if len(percentage_colour) != 6:
        logging.error(
            f" hexadecimal code must contain 6 characters not {len(percentage_colour)}"
        )
        sys.exit(1)

    if not percentage_colour.isalnum():
        logging.error(" hexadecimal code can only contain numbers and letters...")
        sys.exit(1)
    
    transcript_quants = load_transcript_quantification(quants=quantification)
    transcript_list = list(transcript_quants.keys())
	
    if (len(genes_of_interest) > 0) and (len(aliases) > 0):
        for gene in genes_of_interest:
            if gene not in transcript_list:
                logging.error(f" {gene} is not in the list of transcripts...")
                logging.info("Shutting down...")
                sys.exit(1)

    transcript_tpms = list(transcript_quants.values())

    lower = 1 - percentage

    lower_quartile = np.quantile(transcript_tpms, lower)

    counts, bins = np.histogram(
        transcript_tpms,
        bins=50,
    )

    colours = []

    for bin in bins:
        if bin < lower_quartile:
            colours.append("grey")

        else:
            colours.append(f"#{percentage_colour}")

    fig, ax = plt.subplots()

    ax.bar(bins[:-1], counts, width=np.diff(bins), color=colours, edgecolor="none", align="edge")

    # graph customization
    ax.set_xlabel("Transcripts Per Million TPM")  # Set the label for the x-axis
    ax.set_ylabel("Frequency")  # Set the label for the y-axis
    ax.spines["bottom"].set_position(("data", 0))
    ax.set_xlim(
        0,
    )
    ax.spines[["right", "top"]].set_visible(False)

    # Annotate Transcripts of interest on histogram, using aliases if provided
    gene_names = dict(zip(genes_of_interest, genes_of_interest))

    if len(genes_of_interest) == len(aliases):
        gene_names = dict(zip(genes_of_interest, aliases))

    for k, v in gene_names.items():
        if v == "_":
            gene_names[k] = k
    
    i=0
    for gene in genes_of_interest:
        ax.annotate(
            rf"{gene_names[gene]}",
            xy=((transcript_quants[gene], 1)),
            xytext=( (transcript_quants[gene], (counts.max() / 2+i))),
            arrowprops=dict(arrowstyle="-|>", color="black"),
            style="italic",
        )
        i+=2

    # Plot figure legend

    if no_legend is False:
        bottom = mpatches.Patch(color="grey", label=f"Bottom {lower*100}%")
        top = mpatches.Patch(
            color=f"#{percentage_colour}", label=f"Top {percentage*100}%"
        )
        plt.legend(handles=[top, bottom], title="Percentage Expression")

    # Plot and save final figure @ 300 DPI
    plt.tight_layout()
    plt.savefig(f"{graph_out}", dpi=300)

    # Show if toggled on
    if show:
        plt.show()
