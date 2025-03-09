import pandas as pd
import numpy as np
from Bio import SeqIO
import sys

def load_transcripts(transcripts: str) -> list[str]:
    with open(transcripts) as file:
        return [line.strip() for line in file]


def write_transcripts(transcriptome: dict[str,str]) -> None:
    
    for name, seq in transcriptome.items():
        sys.stdout.write(f">{name}\n")
        sys.stdout.write(f"{seq}\n")


def extract_by_id(transcripts: str, transcriptome: str):

    # Load transcript list
    transcripts_list = load_transcripts(transcripts=transcripts)

    # Load transcriptome 
    transcriptome = SeqIO.parse(open(transcriptome), "fasta")

    # Extract sequences for transcripts from transcriptome
    transcripts_of_interest = {record.id : record.seq for record in transcriptome if record.id in transcripts_list}

    # Write extracted sequences to output
    write_transcripts(transcriptome=transcripts_of_interest)

def extract_by_expression(transcripts: str, transcriptome: str, quantification: str, top_percentage: float) -> None: 

    # Load transcript list
    transcripts_list = load_transcripts(transcripts=transcripts)

    # Load transcript quantification
    df = pd.read_csv(quantification, sep="\t")

    # Get expression levels for transcripts in list provided
    transcripts_df = df[(df.Name.isin(transcripts_list))]

    # Find threshold for top X% of transcripts
    threshold = np.quantile(transcripts_df.TPM, (1-top_percentage))

    # Filter for transcripts greater than or equal to threshold
    top_transcripts = transcripts_df[transcripts_df.TPM >= threshold]

    # Get names of transcripts above threshold
    top_transcript_ids = top_transcripts.Name.to_list()

    # Load transcriptome 
    transcriptome = SeqIO.parse(open(transcriptome), "fasta")

    # Extract sequences of transcripts above threshold
    top_transcript_sequences = {record.id : record.seq for record in transcriptome if record.id in top_transcript_ids}

    # Write extracted sequences to output
    write_transcripts(transcriptome=top_transcript_sequences)