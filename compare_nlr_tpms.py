import pandas as pd
from Bio import SeqIO
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches


def load_nlr_transcripts(transcripts: str) -> list[str]:
    with open(transcripts) as file:
        return [line.strip() for line in file]

def visualise_nlr_transcripts(transcripts: str, quantification: str, gene_name: str, graph_out:str, gene_quant: int):
    nlrs = load_nlr_transcripts(transcripts=transcripts)

    df = pd.read_csv(quantification, sep="\t")
    nlrs = df[(df.Name.isin(nlrs))]


    q75 = np.quantile(nlrs.TPM, 0.75)

    counts, bins = np.histogram(nlrs.TPM, bins=150, )

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
    ax.set_ylim(-5, max(counts))
    ax.set_xlim(-2,)
    ax.spines[['right', 'top']].set_visible(False)


    b75 = mpatches.Patch(color='grey', label='Bottom 75%')
    t25 = mpatches.Patch(color='orange', label='Top 25%')

    ax.annotate(fr'{gene_name}', xy=(gene_quant, 1), xytext=(10, 50), textcoords='offset points',
                    arrowprops=dict(arrowstyle='-|>', color='black'), style='italic')


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
    nlr_transcripts = "/home/powellor/Documents/projects/r_gene_expression_david/scripts/nlr_quant_scripts/TA1675_NLR_trancripts.txt"
    bf_transcripts = "/home/powellor/Documents/projects/r_gene_expression_david/scripts/nlr_quant_scripts/TA1675_beta_finger_trancripts.txt"
    quantification = "/home/powellor/Documents/projects/r_gene_expression_david/analysis/TA1675/TA1675_quantification/quant.sf"
    visualise_nlr_transcripts(transcripts=nlr_transcripts, quantification=quantification, gene_name="Yr28", graph_out="Yr28_expression.png", gene_quant=10.195096)
    visualise_nlr_transcripts(transcripts=nlr_transcripts, quantification=quantification, gene_name="Sr46", graph_out="Sr46_TRINITY_DN1954_c0_g1_i12_expression.png", gene_quant=15.449350)
    visualise_nlr_transcripts(transcripts=bf_transcripts, quantification=quantification, gene_name="Lr39", graph_out="Lr39_TRINITY_DN857_c1_g1_i2_expression.png", gene_quant=26.351422)
    visualise_nlr_transcripts(transcripts=bf_transcripts, quantification=quantification, gene_name="WTK4", graph_out="WTK4_TRINITY_DN2349_c0_g2_i1_expression.png", gene_quant=27.023780)
    visualise_nlr_transcripts(transcripts=bf_transcripts, quantification=quantification, gene_name="WTK4", graph_out="WTK4_TRINITY_TRINITY_DN305_c0_g1_i41_expression.png", gene_quant=6.643938)
    """

    """
    TOWWC054
    transcripts = "/home/powellor/Documents/projects/r_gene_expression_david/scripts/nlr_quant_scripts/TOWWC054_NLR_trancripts.txt"
    quantification = "/home/powellor/Documents/projects/r_gene_expression_david/analysis/TOWWC054_quant/TOWWC054_salmon_out/quant.sf"
    visualise_nlr_transcripts(transcripts=transcripts, quantification=quantification, gene_name="SrTA10187", graph_out="SrTA10187_expression.png", gene_quant=13.496783)
    """

    """
    #TA10171
    transcriptome = "/home/powellor/Documents/projects/r_gene_expression_david/analysis/TA10171/trinity_TA10171.Trinity.fasta"
    transcripts = "/home/powellor/Documents/projects/r_gene_expression_david/scripts/nlr_quant_scripts/TA10171_NLR_trancripts.txt"
    quantification = "/home/powellor/Documents/projects/r_gene_expression_david/analysis/TA10171/TA10171_quantification/quant.sf"
    visualise_nlr_transcripts(transcripts=transcripts, quantification=quantification, gene_name="DN2103_c0_g2_i7", graph_out="TA10171_DN2103_c0_g2_i7_expression.png", gene_quant=6.709618)
    visualise_nlr_transcripts(transcripts=transcripts, quantification=quantification, gene_name="Dn2103_c0_g3_i2", graph_out="TA10171_DN2103_c0_g2_i7_expression.png", gene_quant=4.758894)
    """
    
    
    #TOWWC0107
    get_top25(transcriptome="/home/powellor/Documents/projects/r_gene_expression_david/analysis/TOWWC0107/trinity_TOWWC0107.Trinity.fasta", transcripts="/home/powellor/Documents/projects/r_gene_expression_david/scripts/nlr_quant_scripts/TOWWC017_NLR_trancripts.txt", quantification="/home/powellor/Documents/projects/r_gene_expression_david/analysis/TOWWC0107/TOWWC0107_quantification/quant.sf", basename="TOWWC0107_nlrs")
    
    

