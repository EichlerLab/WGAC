#!/usr/bin/env python
                                                                                        
import argparse                                                                         
import sys
                                                                                        
parser = argparse.ArgumentParser()                                                      
parser.add_argument("--szFileOfLegalValues", required = True )                           
parser.add_argument("--n0BasedToken", required = True, type = int )                          
args = parser.parse_args()                                                              


aLegalValues = {}
with open( args.szFileOfLegalValues, "r" ) as fLegal:
    aLines = fLegal.readlines()

    for szLine in aLines:
        szLegalValue = szLine.rstrip()
        aLegalValues[ szLegalValue ] = 1



while True:
    szLine = sys.stdin.readline()
    if ( szLine == "" ):
        break

    aWords = szLine.split()

    if ( len( aWords ) <= args.n0BasedToken ):
        continue

    szWord = aWords[ args.n0BasedToken ]

    #sys.stderr.write( "looking for #" + szWord + "#\n" )


    if ( szWord in aLegalValues ):
        sys.stdout.write( szLine )


