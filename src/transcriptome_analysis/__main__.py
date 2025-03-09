import argparse
import os
import sys
from __init__ import __version__
from visualise_transcript_expression import visualise_transcripts_expression
from extract_transcripts import extract_by_expression, extract_by_id

def is_valid_file(parser, arg):
    if not os.path.exists(arg):
        parser.error(f"Input file ({arg}) not found!")

    else:
        return arg


def main():
    parser = argparse.ArgumentParser(prog="Transcript Analysis Toolkit")

    sub_parsers = parser.add_subparsers(dest="command")

    # Versioning
    parser.add_argument(
        '-v',
        '--version', 
        action='version', 
        version='%(prog)s {version}'.format(version=__version__)
    )

    # sequence extraction related parameters

    extractor_parent_parser = argparse.ArgumentParser(add_help=False)

    extractor_parent_parser.add_argument(
        '-l',
        '--transcript-list',
        metavar='\b', 
        type=lambda x: is_valid_file(parser, x),
        help="List of transcripts",
        required=True,
    )

    extractor_parent_parser.add_argument(
        '-t',
        '--transcriptome',
        metavar='\b', 
        type=lambda x: is_valid_file(parser, x),
        help="Transcriptome",
        required=True,
    )

    # Expression related parameters
    quantification_parent_parser = argparse.ArgumentParser(add_help=False)

    quantification_parent_parser.add_argument(
        "-p",
        "--percentage",
        metavar='\b',
        default=0.25,
        type=float,
        required=True,
        help="Fraction of top percent coloured between 0.0-1.0 (default: 0.25)"
    )

    quantification_parent_parser.add_argument(
        '-q',
        '--quantification',
        metavar='\b', 
        type=lambda x: is_valid_file(parser, x),
        help="Transcript quantification values",
        required=True,
    )

    # Visualiser command and specific parameters

    visualiser_parser = sub_parsers.add_parser("visualise", parents=[quantification_parent_parser])

    visualiser_parser.add_argument(
        '-l',
        '--transcript-list',
        metavar='\b', 
        type=lambda x: is_valid_file(parser, x),
        help="List of transcripts",
        required=True,
    )

    visualiser_parser.add_argument(
        '-o',
        '--output', 
        metavar='\b',
        type=str,
        help="Path to output graph",
        required=True,
    )

    visualiser_parser.add_argument(
        '-i',
        '--transcripts-of-interest', 
        metavar='\b',
        nargs="+",
        type=str,
        help="Names of transcripts to highlight on graph with arrows",
        required=False,
    )

    visualiser_parser.add_argument(
        '-a',
        '--aliases', 
        metavar='\b',
        nargs="+",
        type=str,
        help="Aliases for transcripts of interest. E.g gene names",
        required=False,
    )

    visualiser_parser.add_argument(
        "-c",
        "--colour",
        metavar='\b',
        default="ffa500",
        type=str,
        required=False,
        help="Hexadecimal code to colour the top percentage of transcripts (default: ffa500 [orange])."
    )  

    visualiser_parser.add_argument(
        "-s",
        "--show", 
        action="store_true", 
        help="Show graph upon creation (creates a new window)"
    ) 

    # Extract-by-id command
    sub_parsers.add_parser("extract-by-id", parents=[extractor_parent_parser])

    # Extract-by-threshold command
    sub_parsers.add_parser("extract-by-threshold", parents=[extractor_parent_parser, quantification_parent_parser])

    if len(sys.argv) == 1:
        parser.print_help()
        sys.exit(1)

    args = parser.parse_args()

    if args.command == "visualise":

        visualise_transcripts_expression(
            transcripts=args.transcript_list,
            quantification=args.quantification,
            graph_out=args.output,
            aliases=args.aliases,
            genes_of_interest=args.transcripts_of_interest,
            show=args.show,
            percentage=args.percentage,
            percentage_colour=args.colour
        )
    
    if args.command == "extract-by-id":
        extract_by_id(transcripts=args.transcript_list,
                      transcriptome=args.transcriptome)
    
    if args.command == "extract-by-threshold":

        extract_by_expression(transcripts=args.transcript_list,
                              transcriptome=args.transcriptome,
                              quantification=args.quantification,
                              top_percentage=args.percentage,
                              )

    


if __name__ == "__main__":
    main()
    sys.exit(0)
