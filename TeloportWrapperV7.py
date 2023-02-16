import os
import time
import argparse

#Command line argument set up
parser = argparse.ArgumentParser(
    prog='DenovoTelomereMiner',
    argument_default='try --simple',
    description=(
        'An all-in-one pipeline and anaylsis tool for de-novo Telomeres'
    ),
    epilog=(
        'The -s and -i options are mutally exclusive. '
        'If simple mode is not used options -b is required.'
    ),
)

#command line option for seperated files 
parser.add_argument(
    '-s',
    action='store',
    nargs=2,
    help=(
        'Allows for seperated fastq file names to be inputted on the command line.' 
        'This argument exepcts two file names seperated by a space.'
    ),
)

#command line option for interleaved files
parser.add_argument(
    '-i',
    action='store',
    help=(
        'Allows for the interleaved fastq name to be inputted on the command line.'
        'This argument expects one file name.'
    ),
)

#command line argument for a genome file
parser.add_argument(
    '-g',
    action='store',
    help=(
        'Allows for the assembled genome file name to be inputted on the command line.'
    ),
)

#command line argument for a directory
parser.add_argument(
    '-d',
    action='store',
    help=(
        'The name of the directory where all output files will be placed.'
        'The directory name will default to the name of the assembled genome file'
        'Directory names must be unique'
    ),
)



args = parser.parse_args()

print(args.s)
