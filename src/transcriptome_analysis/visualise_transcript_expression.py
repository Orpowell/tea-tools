import random
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import seaborn as sns
import csv
import logging
import sys
from matplotlib.colors import ListedColormap

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
    percentage: float = 25.0,
    percentage_colour: str = "ffa500",
    no_legend: bool = False,
    svg: bool = False
):
    if (percentage > 100) or (percentage < 0):
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
    
    ax = sns.swarmplot(transcript_tpms, orient="h", color='red', size=3)
    
    cmap = ListedColormap(['grey', f'#{percentage_colour}'])
    
    lower = 100-percentage
    
    for col in ax.collections:
    	x = col.get_offsets()[:,0]
    	perc = np.percentile(x, [lower])
    	col.set_cmap(cmap)
    	col.set_array(np.digitize(x, perc))
    
    sns.violinplot(transcript_tpms, orient="h", inner=None, color='white', linecolor='black', cut=0, density_norm='count')
    ax.set_xlabel("Transcripts Per Million (tpm)")  # Set the label for the x-axis
    ax.spines[["right", "top", "left"]].set_visible(False)
    ax.yaxis.set_ticklabels([])
    ax.yaxis.set_major_formatter(plt.NullFormatter())
    plt.yticks([])

    # Annotate Transcripts of interest on histogram, using aliases if provided
    gene_names = dict(zip(genes_of_interest, genes_of_interest))

    if len(genes_of_interest) == len(aliases):
        gene_names = dict(zip(genes_of_interest, aliases))

    for k, v in gene_names.items():
        if v == "_":
            gene_names[k] = k
    
    x,y = np.array(ax.collections[0].get_offsets()).T
    data_coords = dict(zip(x,y)) 
    
    for gene in genes_of_interest:
    	x = transcript_quants[gene]
    	ax.annotate(
			rf"{gene_names[gene]}",
            xy=((x, data_coords[x])),
    		xytext=((round(x), random.uniform(0.2,0.4))),
	    	arrowprops=dict(arrowstyle="-|>", color="black"),
        	style="italic",
        	)

    # Plot figure legend

    if no_legend is False:
        bottom = mpatches.Patch(color="grey", label=f"Bottom {lower}%")
        top = mpatches.Patch(
            color=f"#{percentage_colour}", label=f"Top {percentage}%"
        )
        plt.legend(handles=[top, bottom], title="Percentage Expression")
    	
    # Plot and save final figure @ 300 DPI
    plt.tight_layout()
    plt.savefig(f"{graph_out}", dpi=300)

    # Show if toggled on
    if show:
        plt.show()
