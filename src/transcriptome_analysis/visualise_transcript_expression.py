import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import csv
import logging
import sys

def load_transcript_list(transcripts: str) -> list[str]:
    
    with open(transcripts) as file:
        return [line.strip() for line in file]

def load_transcript_quantification(quants: str) -> dict[str, float]:

    with open(quants) as file:

        reader = csv.reader(file, delimiter="\t")

        next(reader)

        return {row[0]:float(row[3]) for row in reader }



def visualise_transcripts_expression(transcripts: str, quantification: str, graph_out:str, genes_of_interest=[str], aliases=[str], show=False, percentage: float = 0.25, percentage_colour: str = "ffa500", no_legend: bool=False):

    if (percentage > 1) or (percentage < 0):
        logging.error(f" percentage ({percentage}) must be between 0 and 1.")
        sys.exit(1)
    
    if len(genes_of_interest) != len(aliases):
        logging.error(f" alias ({len(aliases)}) and transcript-of-interest ({len(genes_of_interest)}) lists are not the same length!")
        logging.info("Shutting down...")
        sys.exit(1)
    
    if len(percentage_colour) != 6:
        logging.error(f" hexadecimal code must contain 6 characters not {len(percentage_colour)}")
        sys.exit(1)
    
    if not percentage_colour.isalnum():
        logging.error(" hexadecimal code can only contain numbers and letters...")
        sys.exit(1)
    
    transcript_list = load_transcript_list(transcripts=transcripts)

    for gene in genes_of_interest:
        if gene not in transcript_list:
            logging.error(f" {gene} is not in the list of transcripts...")
            logging.info("Shutting down...")
            sys.exit(1)

    transcript_quants = load_transcript_quantification(quants=quantification)
    
    transcripts_of_interest = {k:v for k,v in transcript_quants.items() if k in transcript_list}

    transcript_tpms = list(transcripts_of_interest.values())

    lower = (1-percentage)

    lower_quartile = np.quantile(transcript_tpms, lower)

    counts, bins = np.histogram(transcript_tpms, bins=150, )

    colours = []

    for bin in bins:
        if bin < lower_quartile:
            colours.append("grey")
        
        else:
            colours.append(f"#{percentage_colour}")

    fig, ax = plt.subplots()

    ax.bar(
        bins[:-1], 
        counts, 
        width=np.diff(bins), 
        color=colours, 
        edgecolor='none'
    ) 

    # graph customization
    ax.set_xlabel('Transcripts Per Million TPM')  # Set the label for the x-axis
    ax.set_ylabel('Frequency')  # Set the label for the y-axis
    ax.spines['bottom'].set_position(('data', 0))
    ax.set_xlim(0,)
    ax.spines[['right', 'top']].set_visible(False)


    # Annotate Transcripts of interest on histogram, using aliases if provided
    gene_names = dict(zip(genes_of_interest, genes_of_interest))

    if len(genes_of_interest) == len(aliases):
        gene_names = dict(zip(genes_of_interest, aliases))
    
    i = 50
    for gene in genes_of_interest:

        ax.annotate(fr'{gene_names[gene]}',
                    xy=(transcript_quants[gene], 1), 
                    xytext=(i, i), 
                    textcoords='offset points',
                    arrowprops=dict(arrowstyle='-|>', color='black'), 
                    style='italic')
        i+=50

    # Plot figure legend
    
    if no_legend is False:
        
        bottom = mpatches.Patch(color='grey', label=f'Bottom {lower*100}%')
        top = mpatches.Patch(color=f'#{percentage_colour}', label=f'Top {percentage*100}%')
        plt.legend(handles=[top, bottom], title="NLR Expression")
    
    # Plot and save final figure @ 300 DPI
    plt.tight_layout()
    plt.savefig(f"{graph_out}", dpi=300)

    # Show if toggled on
    if show:
        plt.show()

        
"""
if __name__ == "__main__":
    

    #visualise_nlr_transcripts(transcripts=nlr_transcripts, quantification=quantification, gene_name="Yr28", graph_out="Yr28_expression.png", gene_quant=10.195096)
    #visualise_nlr_transcripts(transcripts=nlr_transcripts, quantification=quantification, gene_name="Sr46", graph_out="Sr46_TRINITY_DN1954_c0_g1_i12_expression.png", gene_quant=15.449350)

    TOWWC054
    transcripts = "/home/powellor/Documents/projects/r_gene_expression_david/scripts/nlr_quant_scripts/TOWWC054_NLR_trancripts.txt"
    quantification = "/home/powellor/Documents/projects/r_gene_expression_david/analysis/TOWWC054_quant/TOWWC054_salmon_out/quant.sf"
    visualise_nlr_transcripts(transcripts=transcripts, 
                              quantification=quantification, 
                              gene_name="SrTA10187", 
                              graph_out="SrTA10187_expression.png",
                              gene_quant=13.496783)


    
    #TA10171
    visualise_nlr_transcripts(transcripts = "/home/powellor/Documents/projects/r_gene_expression_david/scripts/nlr_quant_scripts/TA10171_NLR_trancripts.txt", 
                              quantification="/home/powellor/Documents/projects/r_gene_expression_david/analysis/TA10171/TA10171_quantification/quant.sf", 
                              graph_out="TA10171_expression.test.png",
                              genes_of_interest=["TRINITY_DN2103_c0_g2_i7","TRINITY_DN2103_c0_g3_i2"],
                              )
    
    
    
    #TOWWC0107
    visualise_nlr_transcripts(transcripts="/home/powellor/Documents/projects/r_gene_expression_david/scripts/nlr_quant_scripts/TOWWC0107_NLR_trancripts.txt", 
                            quantification="/home/powellor/Documents/projects/r_gene_expression_david/analysis/TOWWC0107/TOWWC0107_quantification/quant.sf", 
                            graph_out="TOWWC0107_candidates.png", 
                            genes_of_interest=["TRINITY_DN1105_c0_g1_i7", "TRINITY_DN1105_c0_g1_i6"],
                            )
    """