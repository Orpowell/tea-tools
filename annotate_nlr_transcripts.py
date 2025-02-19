from Bio import SeqIO
import pyhmmer
import pandas as pd
import os
import multiprocessing
import glob
from functools import partial
import shutil
import numpy as np


transcriptome: str = "path"


def transcript_splitter(transcriptome: str, tmp: str) -> None:
    
    records = SeqIO.parse(open(transcriptome), "fasta")

    for record in records:
        with open(f"{tmp}/{record.id}.fasta", "w+") as file:
            file.write(f">{record.id}\n")
            file.write(f"{str(record.seq)}\n")


class Transcript:
    def __init__(self, chunk_path: str, tmp: str) -> None:
        record = next(SeqIO.parse(chunk_path, "fasta"))

        self.transcript_name = record.id
        self.forward_sequence = record.seq
        self.reverse_sequence = record.seq.reverse_complement()
        self.tmp = tmp
        self.forward = f"{self.tmp}/{self.transcript_name}_forward.fasta"
        self.reverse = f"{self.tmp}/{self.transcript_name}_reverse.fasta"

    def translate_ORFs(self):
        with open(self.forward, "w+") as out:
            out.write(">f-0\n")
            out.write(f"{self.forward_sequence.translate()}\n")
            out.write(">f-1\n")
            out.write(f"{self.forward_sequence[1:].translate()}\n")
            out.write(">f-2\n")
            out.write(f"{self.forward_sequence[2:].translate()}\n")

        with open(self.reverse, "w+") as out:
            out.write(">r-0\n")
            out.write(f"{self.reverse_sequence.translate()}\n")
            out.write(">r-1\n")
            out.write(f"{self.reverse_sequence[1:].translate()}\n")
            out.write(">r-2\n")
            out.write(f"{self.reverse_sequence[2:].translate()}\n")

    def run_pyhmmer_hmmsearch(self, input, pfam) -> list[pyhmmer.plan7.Hit]:
        pressed = [f"{pfam}.h3p", f"{pfam}.h3m", f"{pfam}.h3i", f"{pfam}.h3f"]

        if not all(map(os.path.isfile, pressed)):
            pfam_load = pyhmmer.plan7.HMMFile(pfam)
            pyhmmer.hmmer.hmmpress(pfam_load, pfam)

        with pyhmmer.easel.SequenceFile(input, digital=True) as seq_file:
            sequences = list(seq_file)

        with pyhmmer.plan7.HMMFile(pfam) as hmm_file:
            targets = hmm_file.optimized_profiles()

            results = list(
                pyhmmer.hmmsearch(sequences=sequences, queries=targets, cpus=4)
            )

        return results

    def generate_annotation_info(
        self, results: list[pyhmmer.plan7.Hit], strand: bool
    ) -> pd.DataFrame:
        data_array = []

        # Forward strand processing
        if strand:
            for hits in results:
                for hit in hits.reported:
                    for domain in hit.domains.reported:

                        name = domain.alignment.hmm_name.decode()
                        score = domain.i_evalue
                        
                        row = [name, score]

                        data_array.append(row)

        # Reverse strand processing
        else:
            for hits in results:
                for hit in hits.reported:
                    for domain in hit.domains.reported:

                        name = domain.alignment.hmm_name.decode()
                        score = domain.i_evalue

                        row = [name, score]

                        data_array.append(row)

        return pd.DataFrame(data_array)

    def clean_up(self):
        os.remove(self.forward)
        os.remove(self.reverse)

    def run(self, PFAM):
        self.translate_ORFs()

        # Forward
        forward_hits = self.run_pyhmmer_hmmsearch(self.forward, PFAM)
        forward_annotations = self.generate_annotation_info(forward_hits, True)

        # Reverse
        reverse_hits = self.run_pyhmmer_hmmsearch(self.reverse, PFAM)
        reverse_annotations = self.generate_annotation_info(reverse_hits, False)

        # Merge forward and reverse strand annotations
        df = pd.concat([forward_annotations, reverse_annotations], axis=0)

        if len(df) == 0:
            print(f"No annotations in : {self.transcript_name}")

        else:

         df.to_csv(
                f"{self.tmp}/{self.transcript_name}.csv",
                sep="\t",
                index=False,
                header=False,
            )

        self.clean_up()



def analyse_transcript(PFAM: str, tmp: str, chunk_path: str) -> None:
    try:
        transcript = Transcript(chunk_path=chunk_path, tmp=tmp)
        transcript.run(PFAM)
    except ValueError:
        print(f"WARNING! Could not process {chunk_path} - likely Alphabet error...")

def extract_NLRs(tmp: str, file_out: str) -> None:

    transcript_pool = glob.glob(rf"{tmp}/*.csv") 

    with open(f"{file_out}", "w+") as output:
        for transcript in transcript_pool:
            
            name = transcript.split(".")[0].split("/")[-1]
            nbarc_present = False
            lrr_present = False
            
            with open(transcript) as file:
                for line in file:

                    if nbarc_present and lrr_present:
                        break

                    if "NB-ARC" in line:
                        nbarc_present = True

                    if "LRR" in line:
                        lrr_present = True 

            is_nlr = nbarc_present and lrr_present

            if is_nlr:
                output.write(f"{name}\n")


def extract_NLRs(tmp: str, file_out: str) -> None:

    transcript_pool = glob.glob(rf"{tmp}/*.csv") 

    with open(f"{file_out}", "w+") as output:
        for transcript in transcript_pool:
            
            name = transcript.split(".")[0].split("/")[-1]
            nbarc_present = False
            lrr_present = False
            
            with open(transcript) as file:
                for line in file:

                    if nbarc_present and lrr_present:
                        break

                    if "NB-ARC" in line:
                        nbarc_present = True

                    if "LRR" in line:
                        lrr_present = True 

            is_nlr = nbarc_present and lrr_present

            if is_nlr:
                output.write(f"{name}\n")

def extract_BFs(tmp: str, file_out: str) -> None:

    transcript_pool = glob.glob(rf"{tmp}/*.csv") 

    with open(f"{file_out}", "w+") as output:
        for transcript in transcript_pool:
            
            name = transcript.split(".")[0].split("/")[-1]
            BF_present = False
            
            with open(transcript) as file:
                for line in file:

                    if BF_present:
                        break

                    if "finger" in line:
                        BF_present= True

            if BF_present:
                output.write(f"{name}\n")

def annotate_transcriptome(transcriptome: str, outdir: str, cores: int, pfam: str, file_out:str, BF=False) -> None:

    tmp = f"{outdir}/tmp"
    os.mkdir(tmp)

    transcript_splitter(transcriptome=transcriptome, tmp=tmp)

    if cores > multiprocessing.cpu_count():
        print(
            f"Too many cores specified {cores}! {multiprocessing.cpu_count()} cores will be used"
        )
        cores = multiprocessing.cpu_count()

    transcript_pool = glob.glob(rf"{tmp}/*.fasta")

    pool = multiprocessing.get_context("spawn").Pool(processes=cores)
    pool.map(partial(analyse_transcript, pfam, tmp), transcript_pool)

    if BF:
        extract_BFs(tmp=tmp, file_out=file_out)

        try:
            shutil.rmtree(tmp)
        except OSError as e:
            print(f"Error: {e.filename} - {e.strerror}.")

    else:
        extract_NLRs(tmp=tmp, file_out=file_out)

        try:
            shutil.rmtree(tmp)
        except OSError as e:
            print(f"Error: {e.filename} - {e.strerror}.")


if __name__ == "__main__":

    nlr_pfam="/home/powellor/Documents/projects/r_gene_expression_david/scripts/nlr_quant_scripts/NLR.hmm"
    bf_pfam="/home/powellor/Documents/projects/global_wtk/beta_finger_profile_hmm/analysis/chinese_spring_analysis/beta_finger_refinement/refined_beta_finger_hmm/beta_finger_refined.cs.source.hmm"
    
    """
    transcriptome="/home/powellor/Documents/projects/r_gene_expression_david/data/trinity_TOWWC054.Trinity-GG.fasta"
    outdir="/home/powellor/Documents/projects/r_gene_expression_david/scripts/nlr_quant_scripts"
    

    annotate_transcriptome(transcriptome=transcriptome, outdir=outdir, cores=40, pfam=pfam, file_out="TOWWC054_NLR_trancripts.txt")

    #TA10171

    transcriptome="/home/powellor/Documents/projects/r_gene_expression_david/analysis/TA10171/trinity_TA10171.Trinity.fasta"
    outdir="/home/powellor/Documents/projects/r_gene_expression_david/scripts/nlr_quant_scripts"
    pfam="/home/powellor/Documents/projects/r_gene_expression_david/scripts/nlr_quant_scripts/NLR.hmm"

    annotate_transcriptome(transcriptome=transcriptome, outdir=outdir, cores=40, pfam=pfam, file_out="TA10171_NLR_trancripts.txt")
    """
    
    """
    #TA1675
    
    transcriptome="/home/powellor/Documents/projects/r_gene_expression_david/analysis/TA1675/trinity_TA1675.fasta.Trinity.fasta"
    outdir="/home/powellor/Documents/projects/r_gene_expression_david/scripts/nlr_quant_scripts/"

    annotate_transcriptome(transcriptome=transcriptome, outdir=outdir, cores=40, pfam=nlr_pfam, file_out="TA1675_NLR_trancripts.txt")
    annotate_transcriptome(transcriptome=transcriptome, outdir=outdir, cores=40, pfam=bf_pfam, file_out="TA1675_beta_finger_trancripts.txt", BF=True)
    """
    
    #TOWWC0107
    transcriptome="/home/powellor/Documents/projects/r_gene_expression_david/analysis/TOWWC0107/trinity_TOWWC0107.Trinity.fasta"
    outdir="/home/powellor/Documents/projects/r_gene_expression_david/scripts/nlr_quant_scripts"

    annotate_transcriptome(transcriptome=transcriptome, outdir=outdir, cores=40, pfam=nlr_pfam, file_out="TOWWC0107_NLR_trancripts.txt")

