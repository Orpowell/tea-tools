import argparse
import os
import sys
from .__init__ import __version__
from .visualise_transcript_expression import visualise_transcripts_expression

def is_valid_file(parser, arg):
    if not os.path.exists(arg):
        parser.error(f"Input file ({arg}) not found!")

    else:
        return arg


def main():
    parser = argparse.ArgumentParser(prog="Transcript Expression Analysis (T.E.A) Tools")

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
        '-t',
        '--transcriptome',
        metavar='\b', 
        type=lambda x: is_valid_file(parser, x),
        help="Transcriptome",
        required=True
    )

    # Expression related parameters
    quantification_parent_parser = argparse.ArgumentParser(add_help=False)

    quantification_parent_parser.add_argument(
        "-p",
        "--percentage",
        metavar='\b',
        default=25,
        type=float,
        required=False,
        help="Fraction of top percent coloured between 1.0-100.0 (default: 25)"
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
        required=False,
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
        default=[],
        help="Names of transcripts to highlight on graph with arrows",
        required=False,
    )

    visualiser_parser.add_argument(
        '-a',
        '--aliases', 
        metavar='\b',
        nargs="+",
        type=str,
        default=[],
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

    visualiser_parser.add_argument(
        "-n",
        "--no-legend", 
        action="store_true",
        default=False, 
        help="Omit legend from graph"
    )

    if len(sys.argv) == 1:
        parser.print_help()
        sys.exit(1)

    args = parser.parse_args()

    if args.command == "visualise":

        visualise_transcripts_expression(
            quantification=args.quantification,
            graph_out=args.output,
            aliases=args.aliases,
            genes_of_interest=args.transcripts_of_interest,
            show=args.show,
            percentage=args.percentage,
            percentage_colour=args.colour,
            no_legend=args.no_legend
        )
        
    else:
        parser.print_help()
        sys.exit(1)


if __name__ == "__main__":
    main()
    sys.exit(0)
