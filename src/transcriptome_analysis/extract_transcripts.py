import pandas as pd
import numpy as np
from Bio import SeqIO
import sys
import csv

def load_transcripts(transcripts: str) -> list[str]:
    with open(transcripts) as file:
        return [line.strip() for line in file]


def write_transcripts(transcriptome: dict[str, str]) -> None:
    
    for name, seq in transcriptome.items():
        sys.stdout.write(f">{name}\n")
        sys.stdout.write(f"{seq}\n")

def load_transcript_quantification(quants: str) -> dict[str, float]:
    with open(quants) as file:
        reader = csv.reader(file, delimiter="\t")

        return {row[0]: float(row[1]) for row in reader}

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

    # Get expression levels for transcripts in list provided
    transcript_quants = {k:v for k,v in load_transcript_quantification(quants=quantification).items() if k in transcripts_list}

    transcript_tpms = list(transcript_quants.values())

    # Find threshold for top X% of transcripts

    threshold = np.quantile(transcript_tpms, (1-top_percentage))

    # Filter for transcripts greater than or equal to threshold
    top_transcripts = {k: v for k, v in transcript_quants.items() if v >= threshold}

    # Get names of transcripts above threshold
    top_transcript_ids = list(top_transcripts.keys())

    # Load transcriptome 
    transcriptome = SeqIO.parse(open(transcriptome), "fasta")

    # Extract sequences of transcripts above threshold
    top_transcript_sequences = {record.id: record.seq for record in transcriptome if record.id in top_transcript_ids}

    # Write extracted sequences to output
    write_transcripts(transcriptome=top_transcript_sequences)