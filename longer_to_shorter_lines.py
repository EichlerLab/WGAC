#!/usr/bin/env python


# module load miniconda/4.5.12

from Bio import SeqIO

import argparse

parser = argparse.ArgumentParser()
parser.add_argument("--szInputFile", required = True )
parser.add_argument("--szOutputFile", required = True )
args = parser.parse_args()


szInput = args.szInputFile
szOutput = args.szOutputFile


with open( szInput, "r" ) as fInput, open( szOutput, "w" ) as fOutput:

    
    iterInput = SeqIO.parse( fInput, "fasta")
    iterOutput = ( record for record in iterInput )

    SeqIO.write( iterOutput, fOutput, "fasta")


