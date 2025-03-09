import argparse
import os
import sys
from __init__ import __version__
from visualise_transcript_expression import visualise_transcripts_expression

def is_valid_file(parser, arg):
    if not os.path.exists(arg):
        parser.error(f"Input file ({arg}) not found!")

    else:
        return arg


def main():
    parser = argparse.ArgumentParser(prog="Transcript Analysis Toolkit")

    sub_parsers = parser.add_subparsers(dest="command")

    visualiser = sub_parsers.add_parser("visualise")


    #############
    # Visualiser #
    #############
    visualiser.add_argument(
        '-v',
        '--version', 
        action='version', 
        version='%(prog)s {version}'.format(version=__version__)
    )

    visualiser.add_argument(
        '-l',
        '--transcripts-list',
        metavar='\b', 
        type=lambda x: is_valid_file(parser, x),
        help="List of transcripts",
        required=True,
    )

    visualiser.add_argument(
        '-q',
        '--quantification',
        metavar='\b', 
        type=lambda x: is_valid_file(parser, x),
        help="Transcript quantification values",
        required=True,
    )

    visualiser.add_argument(
        '-o',
        '--output', 
        metavar='\b',
        type=str,
        help="Path to output graph",
        required=True,
    )

    visualiser.add_argument(
        '-i',
        '--transcripts-of-interest', 
        metavar='\b',
        nargs="+",
        type=str,
        help="Names of transcripts to highlight on graph with arrows",
        required=False,
    )

    visualiser.add_argument(
        '-a',
        '--aliases', 
        metavar='\b',
        nargs="+",
        type=str,
        help="Aliases for transcripts of interest. E.g gene names",
        required=False,
    )

    visualiser.add_argument(
        "-p",
        "--percentage",
        metavar='\b',
        default=0.25,
        type=float,
        required=False,
        help="Fraction of top percent coloured between 0.0-1.0 (default: 0.25)"
    )

    visualiser.add_argument(
        "-pc",
        "--percentage-colour",
        metavar='\b',
        default="ffa500",
        type=str,
        required=False,
        help="Hexadecimal code to colour the top percentage of transcripts (default: ffa500 [orange])."
    )  

    visualiser.add_argument(
        "-s",
        "--show", 
        action="store_true", 
        help="Show graph upon creation (creates a new window)"
    ) 
        
    if len(sys.argv) == 1:
        parser.print_help()
        sys.exit(1)

    args = parser.parse_args()

    if args.command == "visualise":

        visualise_transcripts_expression(
            transcripts=args.transcripts_list,
            quantification=args.quantification,
            graph_out=args.output,
            aliases=args.aliases,
            genes_of_interest=args.transcripts_of_interest,
            show=args.show,
            percentage=args.percentage,
            percentage_colour=args.percentage_colour
        )


if __name__ == "__main__":
    main()
    sys.exit(0)
