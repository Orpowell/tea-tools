# Transcriptome Expression Analysis (T.E.A) Tools

Tea-tools is a toolkit for investigating transcripts expression in transcriptomes quantified using salmon.  

## Installation

Tea-tools has only been tested using Python 3.10.14. Therefore, we recommend installing tea-tools in an environment using conda or venv.

### Download tea-tools

    git clone https://github.com/orpowell/tea-tools.git
    cd tea-tools

### Installation with conda (recommended)

    conda create -n tea-tools python=3.10.14
    conda activate tea-tools
    pip install .

# Preparing a Transcript Quantification Input for Tea-tools

Several of Tea-tools commands require inforamtion about transcript quantification that can be generated from a number of tools. For tea-tools transcript quantification data should be provided in the format of a 2 column headerless TSV file. The first column is the transcript name. The second column is the transcript TPM value. 

This data can be easily extracted from the output of quantification tools using a mixture of AWK and the bash command tail. An example using Salmon is shown below:

    salmon index -t transcriptome.fna -i transcriptome_index

    salmon quant -i transcriptome_index \
            -l A \
            -1 reads_1.fastq.gz \
            -2 reads_2.fastq.gz \
            -o transcript_quantification

    # The transcript quantification file is found in transcript_quantification/quant.sf

    awk '{print $1"\t"$4}' transcript_quantification/quant.sf | tail -n +2 > transcript-quantification.tsv


# Running Tea-tools

Three commands are currently available using tea-tools: visualise, extract-by-id, extract-by-threshold. All commands can be run using the following command:

    tea-tools <command> ...

# Visualise

The visualise command generates a histogram of transcript expression for set of specified transcripts and highlights the top x% most expressed transcripts. By default, the top 25% of transcripts are highlighted. Specific transcripts can be indicated with an arrow and labelled if desired (See: parameters table below).

Run the visualise command as follows:

    tea-tools visualise \
        --quantification transcript-quantification.tsv \
        --transcripts-of-interest TRANSCRIPT_1 TRANSCRIPT_2 \
        --aliases GENE_1 _ \
        --percentage 0.25 \
        -show \
        --output transcript_expression.png

## Parameters

|Parameter|Description|
|---|---|
|--quantification| Path to transcript quantification data (See: Preparing a Transcript Quantification Input for Tea-tools).|
|--transcripts-of-interest| Individual transcript(s) to be indicated on the histogram (Optional). |
|--aliases | Alternative names for transcripts specified in --transcripts-of-interest (Optional). Aliases should be listed in the same order as the transcripts are listed. To retain the original transcripts name specify "_" instead of an alias.|
|--output| Path where to save PNG file of graph. |
|--percentage| Top percentage most expressed transcripts to colour between 1.0-100.0 (Deafult: 25, I.e top 25% most expressed transcripts). |
|--colour| Hexadecimal code to colour most expressed transcripts (Default: ffa500 [Orange]).|
|--show| Immediately display the generated graph in a pop-up window. |
|--no-legend| Omit the legend from the generated graph. |


# Citing tea-tools

TBA

# Contact

If you wish to discuss tea-tools, raise issues or collaborate please contact me [here](mailto:mail@oliverpowell.com)
