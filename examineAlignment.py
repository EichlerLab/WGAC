#!/usr/bin/env python

# WGAC
# this script will examine one of the alignments in data/align_both and give various 
# measures of the quality of that alignment
# David Gordon, Dec 2022


import argparse


parser = argparse.ArgumentParser()
parser.add_argument("--szAlignment", required = True )
args = parser.parse_args()


with open( args.szAlignment, "r" ) as fAlignment:

    aLines = fAlignment.readlines()


# looks like this:

# FAST ALIGN   chr1 (223987267 to 223992675) vs chr4 (118596805 to 118591686)  line 0
# Global alignment with 5434 spaces. -f -40 -g -1                              line 1
#                                                                               line 2
#      223987270 223987280 223987290 223987300 223987310 223987320             line 3
# chr1      CCTCAGAGCACCCCCAGGTATTCCAACCTAATCCTGGTGCCCCGCCTCTCACCACCCTTC       line 4
# 1         ||||||||||||||||||||||||||||||||||||||*|||||||||||||||||||||       line 5
# chr4      CCTCAGAGCACCCCCAGGTATTCCAACCTAATCCTGGTACCCCGCCTCTCACCACCCTTC       line 6
#        118596800 118596790 118596780 118596770 118596760 118596750           line 7
#                                                                              line 8   
#      223987330  223987340 223987350 223987360 223987370 223987380            line 9
# chr1      TTCCTG-TTTAACCTCAACCCCTACACAAAGCCTGGGCCACTTAATGTGGCATCAAACAG       line 10
# 61        |||||| |||||||||||*|||||||||||||||||||||||||||||||||||||||||
# chr4      TTCCTGCTTTAACCTCAATCCCTACACAAAGCCTGGGCCACTTAATGTGGCATCAAACAG
#        118596740 118596730 118596720 118596710 118596700 118596690


szTopLineBases =""
nLine = 4
while nLine < ( len( aLines ) - 1 ):
    szLine = aLines[ nLine ].rstrip()

    # get rid of leading labels which takes up 10 spaces

    szTopLineBases += szLine[10:]
    # set up for next round
    nLine += 6
    
szBotLineBases = ""
nLine = 6
while nLine < ( len( aLines ) - 1 ):
    szLine = aLines[ nLine ].rstrip()

    # get rid of leading labels which takes up 10 spaces

    szBotLineBases += szLine[10:]
    # set up for next round
    nLine += 6

assert len( szTopLineBases ) == len( szBotLineBases )

nMatches = 0
nMismatches = 0
nGaps = 0
bInGap = False
nLargestGap = 0
nCurrentGapSize = 0

for n in range( 0, len( szTopLineBases ) ):
    cTopBase = szTopLineBases[n]
    cBotBase = szBotLineBases[n]


    bGapHere = False
    
    if (cTopBase == "-" or cBotBase == "-" ):
        nGaps += 1
        bGapHere = True
    elif( cTopBase == cBotBase ):
        nMatches += 1
    else:
        nMismatches += 1


    if bInGap:
        if ( bGapHere ):
            nCurrentGapSize += 1
        else:
            # gap ended
            bInGap = False
            if ( nCurrentGapSize > nLargestGap ):
                nLargestGap = nCurrentGapSize
    else:
        if ( bGapHere ):
            # gap started
            bInGap = True
            nCurrentGapSize = 1
        # otherwise just continue no gap

print( "aligned length: {:d} matches: {:d} mismatches: {:d} largest gap: {:d} gaps: {:d} identity ignoring gaps: {:.1f} % identity counting gaps: {:.1f} % of alignment that is gaps: {:.1f} %".format( len( szTopLineBases ), nMatches, nMismatches, nLargestGap, nGaps, nMatches * 100.0 / ( nMatches + nMismatches ), nMatches * 100.0 / len( szTopLineBases ), nGaps * 100.0 / len( szTopLineBases ) ) )


                                                                                          
