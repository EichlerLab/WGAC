#!/usr/bin/env python


import subprocess
import sys
import argparse

parser = argparse.ArgumentParser()
parser.add_argument("--szSgeAlignBatch", required = True )
parser.add_argument("--szSgeAlignBatch110k", required = True )
parser.add_argument("--global_align_110k_flags_directory", required = True )
parser.add_argument("--fMaximumPerCentAlignmentFailure", required = True )
args = parser.parse_args()                                                                          




with open( args.szSgeAlignBatch, "r" ) as fSgeAlignBatch:
    n = 0
    while True:
        szLine = fSgeAlignBatch.readline()
        if ( szLine == "" ):
            break
        n += 1

nTotalAlignmentsToAttempt = n

with open( args.szSgeAlignBatch110k, "r" ) as fSgeAlignBatch110k:
    n = 0
    while True:
        szLine = fSgeAlignBatch110k.readline()
        if ( szLine == "" ):
            break
        n += 1

n110kAlignmentsToAttempt = n

szCommand = "ls " + args.global_align_110k_flags_directory + " | wc -l"
print( "about to execute: " + szCommand )
szOutput = subprocess.check_output( szCommand, shell = True )

n110kAlignmentsSucceeded = int( szOutput.rstrip() )

nFailed = n110kAlignmentsToAttempt - n110kAlignmentsSucceeded

fPerCentFailed = float( nFailed ) * 100.0 / float( nTotalAlignmentsToAttempt )

print( "alignments attempted: {:n}".format( nTotalAlignmentsToAttempt ) )
print( "alignments failed: {:n} or {:.6f} %".format( nFailed, fPerCentFailed  ) )


if ( fPerCentFailed < float( args.fMaximumPerCentAlignmentFailure ) ):
    sys.stderr.write("acceptable number of failures so pipeline is continuing\n")
    sys.exit( 0 )
else:
    sys.stderr.write("too many alignments failed so pipeline is stopping\n")
    sys.exit( 1 )
