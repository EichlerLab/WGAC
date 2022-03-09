#!/usr/bin/env python


import argparse
parser = argparse.ArgumentParser()
parser.add_argument("--szOffsetTable", required = True )
parser.add_argument("--xm_al", required = True )
parser.add_argument("--szOutputxm_al", required = True )
args = parser.parse_args()
                                                                          
dictOffset = {}

nLengthOfAllUnplacedContigs = 0

with open( args.szOffsetTable, "r" ) as fOffsetTable:

    for szLine in fOffsetTable.readlines():

        aWords = szLine.split()
        # looks like:
        # Super_Scaffold_86       5696367 7895226
        # 0                         1       2
        # name                    length   offset

        if ( len( aWords ) != 3 ):
            print "this line didn't have 3 tokens but should: " + szLine
            continue

        nLengthOfContig = int( aWords[1] )
        nLengthOfAllUnplacedContigs += nLengthOfContig

        dictOffset[ aWords[0] ] = int( aWords[2] )

with open( args.xm_al, "r" ) as fInput, open( args.szOutputxm_al, "w" ) as fOutput:
    while True:
        szLine = fInput.readline()
        if ( szLine == "" ):
            break


        if ( szLine.startswith( "QNAME" ) ):
            fOutput.write( szLine )
            continue


        aWords = szLine.split()
        # looks like:
        # QNAME   QB      QE      QLEN    SNAME   SB      SE      SLEN    base_S  per_sim 1
        # 000001F_75230254_qpds_75173938_75177834_scaf    1       3897    3897    000949F_119402_qpd_scaf 23362   19432   120103  1007    0.915590863952334
        # 0       1       2        3        4     5       6        7


        szContig1 = aWords[0]
        szContig2 = aWords[4]

        if ( szContig1 in dictOffset ):
            nOffset = dictOffset[ szContig1 ]
            nStart1 = int( aWords[1] ) + nOffset
            nEnd1 = int( aWords[2] ) + nOffset

            aWords[0] = "UNK"
            aWords[1] = str( nStart1 )
            aWords[2] = str( nEnd1 )
            aWords[3] = str( nLengthOfAllUnplacedContigs )

            
        if ( szContig2 in dictOffset ):
            nOffset = dictOffset[ szContig2 ]
            nStart2 = int( aWords[5] ) + nOffset
            nEnd2 = int( aWords[6] ) + nOffset
            
            aWords[4] = "UNK"
            aWords[5] = str( nStart2 )
            aWords[6] = str( nEnd2 )
            aWords[7] = str( nLengthOfAllUnplacedContigs )


        szLineToWrite = "\t".join( aWords )
        fOutput.write( szLineToWrite + "\n" )
        
            
  
