import pandas as pd
from Bio import SeqIO
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import csv

def load_nlr_transcripts(transcripts: str) -> list[str]:
    with open(transcripts) as file:
        return [line.strip() for line in file]

def load_transcript_quantification(quants: str) -> dict[str, float]:

    with open(quants) as file:

        reader = csv.reader(file, delimiter="\t")

        next(reader)

        return {row[0]:float(row[3]) for row in reader }


def visualise_nlr_transcripts(transcripts: str, quantification: str, graph_out:str, genes_of_interest=[str], aliases=[str]):
    nlrs = load_nlr_transcripts(transcripts=transcripts)

    transcript_quants = load_transcript_quantification(quants=quantification)
    
    nlrs = {k:v for k,v in transcript_quants.items() if k in nlrs}

    nlr_tpm = list(nlrs.values())

    q75 = np.quantile(nlr_tpm, 0.75)

    counts, bins = np.histogram(nlr_tpm, bins=150, )

    colours = []

    for bin in bins:
        if bin < q75:
            colours.append("grey")
        
        else:
            colours.append("orange")

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
    ax.spines['bottom'].set_position(('data', -5))
    #ax.set_ylim(-5, max(counts))
    ax.set_xlim(-2,)
    ax.spines[['right', 'top']].set_visible(False)


    b75 = mpatches.Patch(color='grey', label='Bottom 75%')
    t25 = mpatches.Patch(color='orange', label='Top 25%')

    gene_names = dict(zip(genes_of_interest, genes_of_interest))

    if len(genes_of_interest) == len(aliases):
        gene_names = dict(zip(genes_of_interest, aliases))
    

    i = 50
    for gene in genes_of_interest:


        
        
        ax.annotate(fr'{gene_names[gene]}', xy=(transcript_quants[gene], 1), xytext=(i, i), textcoords='offset points',
                        arrowprops=dict(arrowstyle='-|>', color='black'), style='italic')
        i+=50

    plt.legend(handles=[t25, b75], title="NLR Expression")
    plt.tight_layout()
    plt.savefig(f"{graph_out}", dpi=300)
    plt.show()


      # Extract top25
def get_top25(transcripts: str, transcriptome: str, quantification: str, basename: str) -> None: 
    nlrs = load_nlr_transcripts(transcripts=transcripts)
    records = SeqIO.parse(open(transcriptome), "fasta")

    nlr_transcripts = {record.id : record.seq for record in records if record.id in nlrs}

    with open(f"{basename}_transcripts.fna", "w+") as out:
        for name, seq in nlr_transcripts.items():
            out.write(f">{name}\n")
            out.write(f"{seq}\n")

    df = pd.read_csv(quantification, sep="\t")
    nlrs_df = df[(df.Name.isin(nlrs))]

    print(nlrs_df.describe())

    q75 = np.quantile(nlrs_df.TPM, 0.75)

    print(q75)

    nlrs_top_25 = nlrs_df[nlrs_df.TPM >= q75]

    top25nlrs=nlrs_top_25.Name.to_list()

    records = SeqIO.parse(open(transcriptome), "fasta")
    nlrs_top_25_transcripts = {record.id : record.seq for record in records if record.id in top25nlrs}

    print(nlrs_top_25_transcripts)

    with open(f"{basename}_transcripts.top25.fna", "w+") as out:
        for name, seq in nlrs_top_25_transcripts.items():
            out.write(f">{name}\n")
            out.write(f"{seq}\n")
  
        

if __name__ == "__main__":
    
    
    """
    #visualise_nlr_transcripts(transcripts=nlr_transcripts, quantification=quantification, gene_name="Yr28", graph_out="Yr28_expression.png", gene_quant=10.195096)
    #visualise_nlr_transcripts(transcripts=nlr_transcripts, quantification=quantification, gene_name="Sr46", graph_out="Sr46_TRINITY_DN1954_c0_g1_i12_expression.png", gene_quant=15.449350)

    visualise_nlr_transcripts(transcripts="/home/powellor/Documents/projects/r_gene_expression_david/scripts/nlr_quant_scripts/TA1675_beta_finger_trancripts.txt", 
                              quantification="/home/powellor/Documents/projects/r_gene_expression_david/analysis/TA1675/TA1675_quantification/quant.sf", 
                              graph_out="TA1675_cloned_r_genes.png", 
                              genes_of_interest=["TRINITY_DN857_c1_g1_i2", "TRINITY_DN2349_c0_g2_i1", "TRINITY_DN305_c0_g1_i41"],
                              aliases=["Lr39", "WTK4_1", "WTK4_2"]
                              )     



    # Cloned KFPs in TA1675 
    visualise_nlr_transcripts(transcripts="/home/powellor/Documents/projects/r_gene_expression_david/scripts/nlr_quant_scripts/TA1675_beta_finger_trancripts.txt", 
                              quantification="/home/powellor/Documents/projects/r_gene_expression_david/analysis/TA1675/TA1675_quantification/quant.sf", 
                              graph_out="TA1675_cloned_r_genes.png", 
                              genes_of_interest=["TRINITY_DN857_c1_g1_i2", "TRINITY_DN2349_c0_g2_i1", "TRINITY_DN305_c0_g1_i41"],
                              aliases=["Lr39", "WTK4_1", "WTK4_2"]
                              )    


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