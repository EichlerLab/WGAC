#!/bin/env python

import argparse
 
 
import pandas as pd


parser = argparse.ArgumentParser()
parser.add_argument("--szSummaryInputFile", required = True )
parser.add_argument("--szExcelOutputFile",  required = True )
args = parser.parse_args()

aLines = []
with open( args.szSummaryInputFile, 'r' ) as fInput:
    for szLine in fInput.readlines():
        aLines.append( szLine.rstrip() )
    


# values of nState
nInitial = 1
nBp = 2
nNonredundantLoci = 3
# the following are used only to prevent the state staying at nNonredundantLoci
# which will cause a few values to be overwritten
nNonredundantLociInterChromosomal = 4
nNonredundantLociIntraChromosomal = 5


nState = nInitial

szPairwiseAlignments1kb90Per = "-666"
szPairwiseAlignmentsJustChr1kb90Per = "-666"
szPairwiseAlignmentsJustChrInter1kb90Per = "-666"
szPairwiseAlignmentsJustChrIntra1kb90Per = "-666"
szPairwiseAlignments10kb95Per = "-666"
szPairwiseAlignmentsJustChr10kb95Per = "-666"
szPairwiseAlignmentsJustChrInter10kb95Per = "-666"
szPairwiseAlignmentsJustChrIntra10kb95Per = "-666"
#
szBpNonredundant1kb90Per = "-666"
szBpRedundant1kb90Per = "-666"
szBpNonredundantJustChr1kb90Per = "-666"
szBpRedundantJustChr1kb90Per = "-666"
szBpRedundantJustChrInter1kb90Per = "-666"
szBpRedundantJustChrIntra1kb90Per = "-666"
#
szBpNonredundant10kb95Per = "-666"
szBpRedundant10kb95Per = "-666"
szBpNonredundantJustChr10kb95Per = "-666"
szBpRedundantJustChr10kb95Per = "-666"
szBpRedundantJustChrInter10kb95Per = "-666"
szBpRedundantJustChrIntra10kb95Per = "-666"
#
szNonredundantLoci1kb90Per = "-666"
szNonredundantLociJustChr1kb90Per = "-666"
szNonredundantLoci10kb95Per = "-666"
szNonredundantLociJustChr10kb95Per = "-666"


for nLine in range( 0, len( aLines ) ):
    szLine = aLines[nLine]

    if ( szLine.startswith( "temp1kb_90percent.tab" ) and ( aLines[nLine+1] ).startswith( "# of pairwise alignments:" ) ):
        szPairwiseAlignments1kb90Per = aLines[nLine+2]

    if ( szLine == "temp1kb_90percent_just_chr.tab" and aLines[nLine+1] == "# of pairwise alignments:" ):
        szPairwiseAlignmentsJustChr1kb90Per = aLines[nLine + 2]

        if ( aLines[nLine + 3 ] == "  inter:" ):
            szPairwiseAlignmentsJustChrInter1kb90Per = aLines[ nLine + 4 ]

        if ( aLines[nLine + 5 ] == "  intra:" ):
            szPairwiseAlignmentsJustChrIntra1kb90Per = aLines[ nLine + 6 ]

        # # debugging
        # if ( szLine == "temp10kb_95percent.tab" ):
        #     print( "almost found, next line = " + aLines[ nLine + 1 ] )
        # # end debugging
        

    if ( szLine == "temp10kb_95percent.tab" and aLines[ nLine + 1 ] == "# of pairwise alignments:" ):
        szPairwiseAlignments10kb95Per = aLines[ nLine + 2 ]


    if ( szLine == "temp10kb_95percent_just_chr.tab" and aLines[ nLine + 1 ] == "# of pairwise alignments:" ):
        szPairwiseAlignmentsJustChr10kb95Per = aLines[ nLine + 2 ]

        if ( aLines[ nLine + 3 ] == "  inter:" ):
            szPairwiseAlignmentsJustChrInter10kb95Per = aLines[ nLine + 4 ]
        
        if ( aLines[ nLine + 5 ] == "  intra:" ):
            szPairwiseAlignmentsJustChrIntra10kb95Per = aLines[ nLine + 6 ]

    if ( szLine == "now bp" ):
        nState = nBp

    if ( nState == nBp and szLine == "temp1kb_90percent.tab" and aLines[ nLine + 1 ] == "# of bases (nonredundant, redundant):" ):
        szBpNonredundant1kb90Per = aLines[ nLine + 2 ]
        szBpRedundant1kb90Per    = aLines[ nLine + 3 ]

    if ( nState == nBp and szLine == "temp1kb_90percent_just_chr.tab" and aLines[ nLine + 1 ] == "# of bases (nonredundant, redundant):" ):
        szBpNonredundantJustChr1kb90Per = aLines[ nLine + 2 ]
        szBpRedundantJustChr1kb90Per    = aLines[ nLine + 3 ]

        if ( aLines[ nLine + 4 ] == "  inter (nonredundant, redundant):" ):
            szBpRedundantJustChrInter1kb90Per = aLines[ nLine + 6 ]

    
        if ( aLines[ nLine + 7 ] == "  intra (nonredundant, redundant):" ):
            szBpRedundantJustChrIntra1kb90Per = aLines[ nLine + 9 ]


    if ( nState == nBp and szLine == "temp10kb_95percent.tab" and aLines[ nLine + 1 ] == "# of bases (nonredundant, redundant):" ):
        szBpNonredundant10kb95Per = aLines[ nLine + 2 ]
        szBpRedundant10kb95Per    = aLines[ nLine + 3 ]

    if ( nState == nBp and szLine == "temp10kb_95percent_just_chr.tab" and aLines[ nLine + 1 ] == "# of bases (nonredundant, redundant):" ):
        szBpNonredundantJustChr10kb95Per = aLines[ nLine + 2 ]
        szBpRedundantJustChr10kb95Per    = aLines[ nLine + 3 ]

        if ( aLines[ nLine + 4 ] == "  inter (nonredundant, redundant):" ):
            szBpRedundantJustChrInter10kb95Per = aLines[ nLine + 6 ]

        if ( aLines[ nLine + 7 ] == "  intra (nonredundant, redundant):" ):
            szBpRedundantJustChrIntra10kb95Per = aLines[ nLine + 9 ]

    if ( szLine == "now nonredundant loci" ):
        nState = nNonredundantLoci

    if ( nState == nNonredundantLoci and szLine == "temp1kb_90percent.tab" ):
        szNonredundantLoci1kb90Per = aLines[ nLine + 1 ]
         
    if ( nState == nNonredundantLoci and szLine == "temp1kb_90percent_just_chr.tab" ):
        szNonredundantLociJustChr1kb90Per = aLines[ nLine + 1 ]

    if ( nState == nNonredundantLoci and szLine == "temp10kb_95percent.tab" ):
        szNonredundantLoci10kb95Per = aLines[ nLine + 1 ]

    if ( nState == nNonredundantLoci and szLine == "temp10kb_95percent_just_chr.tab" ):
        szNonredundantLociJustChr10kb95Per = aLines[ nLine + 1 ]

    if ( szLine == "inter-chromosomal" ):
        nState = nNonredundantLociInterChromosomal

    if ( szLine == "intra-chromosomal" ):
        nState = nNonredundantLociIntraChromosomal

             
            
print( "szPairwiseAlignments1kb90Per = " + szPairwiseAlignments1kb90Per )
print( "szPairwiseAlignmentsJustChr1kb90Per = " + szPairwiseAlignmentsJustChr1kb90Per )
print( "szPairwiseAlignmentsJustChrInter1kb90Per = " + szPairwiseAlignmentsJustChrInter1kb90Per )
print( "szPairwiseAlignmentsJustChrIntra1kb90Per = " + szPairwiseAlignmentsJustChrIntra1kb90Per )
print( "szPairwiseAlignments10kb95Per = " + szPairwiseAlignments10kb95Per )
print( "szPairwiseAlignmentsJustChr10kb95Per = " + szPairwiseAlignmentsJustChr10kb95Per )
print( "szPairwiseAlignmentsJustChrInter10kb95Per = " + szPairwiseAlignmentsJustChrInter10kb95Per )
print( "szPairwiseAlignmentsJustChrIntra10kb95Per = " + szPairwiseAlignmentsJustChrIntra10kb95Per )
#
print( "szBpNonredundant1kb90Per = " + szBpNonredundant1kb90Per )
print( "szBpRedundant1kb90Per = " + szBpRedundant1kb90Per )
print( "szBpNonredundantJustChr1kb90Per = " + szBpNonredundantJustChr1kb90Per )
print( "szBpRedundantJustChr1kb90Per = " + szBpRedundantJustChr1kb90Per )
print( "szBpRedundantJustChrInter1kb90Per = " + szBpRedundantJustChrInter1kb90Per )
print( "szBpRedundantJustChrIntra1kb90Per = " + szBpRedundantJustChrIntra1kb90Per )
#
print( "szBpNonredundant10kb95Per = " + szBpNonredundant10kb95Per )
print( "szBpRedundant10kb95Per = " + szBpRedundant10kb95Per )
print( "szBpNonredundantJustChr10kb95Per = " + szBpNonredundantJustChr10kb95Per )
print( "szBpRedundantJustChr10kb95Per = " + szBpRedundantJustChr10kb95Per )
print( "szBpRedundantJustChrInter10kb95Per = " + szBpRedundantJustChrInter10kb95Per )
print( "szBpRedundantJustChrIntra10kb95Per = " + szBpRedundantJustChrIntra10kb95Per )
#
print( "szNonredundantLoci1kb90Per = " + szNonredundantLoci1kb90Per )
print( "szNonredundantLociJustChr1kb90Per = " + szNonredundantLociJustChr1kb90Per )
print( "szNonredundantLoci10kb95Per = " + szNonredundantLoci10kb95Per )
print( "szNonredundantLociJustChr10kb95Per = " + szNonredundantLociJustChr10kb95Per )


df1 = pd.DataFrame([ 
['Counts', 'Base Pairs', 'Counts', 'Base Pairs' ],
['','','',''], 

[szPairwiseAlignmentsJustChr1kb90Per, szBpRedundantJustChr1kb90Per, szPairwiseAlignments1kb90Per, szBpRedundant1kb90Per], 
[szPairwiseAlignmentsJustChrInter1kb90Per, szBpRedundantJustChrInter1kb90Per, 'NA', 'NA'], 
[szPairwiseAlignmentsJustChrIntra1kb90Per, szBpRedundantJustChrIntra1kb90Per, 'NA', 'NA'], 
[szNonredundantLociJustChr1kb90Per, szBpNonredundantJustChr1kb90Per, szNonredundantLoci1kb90Per, szBpNonredundant1kb90Per], 

['','','',''], 


[szPairwiseAlignmentsJustChr10kb95Per, szBpRedundantJustChr10kb95Per, szPairwiseAlignments10kb95Per, szBpRedundant10kb95Per], 
[szPairwiseAlignmentsJustChrInter10kb95Per, szBpRedundantJustChrInter10kb95Per, 'NA', 'NA'], 
[szPairwiseAlignmentsJustChrIntra10kb95Per, szBpRedundantJustChrIntra10kb95Per, 'NA', 'NA'], 
[szNonredundantLociJustChr10kb95Per, szBpNonredundantJustChr10kb95Per, szNonredundantLoci10kb95Per, szBpNonredundant10kb95Per] 
], 
                   index=['', '1kb length, 90% similarity', 'Pairwise Alignments', 'Interchromosome Pair. Align.', 'Intrachromosomal Pair. Align.', 'Nonredundant Loci', '10kb length, 95% similarity:', 'Pairwise Alignments', 'Interchromosome Pair. Align.', 'Intrachromosomal Pair. Align.', 'Nonredundant Loci'],
                   columns=['Just Chromosomes', '', 'Chromosomes + UNK', ''])

df1.to_excel( args.szExcelOutputFile )

print( "see: " + args.szExcelOutputFile )
