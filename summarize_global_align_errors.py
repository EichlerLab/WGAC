#!/usr/bin/env python

import subprocess

with open( "summarize_global_align_errors.txt", "w" ) as fSummarizeOutput:

    szCommand = "ls -tr .snakemake/log/* | tail -1"
    print "about to execute: " + szCommand
    szOutput = subprocess.check_output( szCommand, shell = True )

    szLogFile = szOutput.rstrip()

    szCommand = "grep '^Error executing' " + szLogFile + " | awk  '{print $10}' | sed 's/,//' "
    print "about to execute: " + szCommand
    szOutput = subprocess.check_output( szCommand, shell = True )

    aJobIDs = szOutput.splitlines()

    nBadMpos = 0
    nNoBadMpos = 0


    for szJobID in aJobIDs:
        szFileWithAsterisks = "log/*.e" + szJobID
        szCommand = "ls " + szFileWithAsterisks
        print "about to execute: " + szCommand
        szOutput = subprocess.check_output( szCommand, shell = True )

        szFile = szOutput.rstrip()

        fSummarizeOutput.write( "log file: " + szFile + "\n" )

        bFoundMkdir = False
        nLinesAfterCommand = 0
        szErrorLines = ""

        with open( szFile, "r" ) as fLog:
            for szLine in fLog.readlines():
                if ( szLine.startswith( "mkdir " ) ):
                    # looks like:
                    # mkdir -p /net/eichler/vol26/projects/wssd_wgac/nobackups/dgordon/wgac/panTro6/both_tmp/both0113351 && cd /net/eichler/vol26/projects/wssd_wgac/nobackups/dgordon/wgac/panTro6/both_tmp/both0113351 && /net/eichler/vol26/projects/wssd_wgac/nobackups/dgordon/wgac/panTro6/align_fast3.pl -i /net/eichler/vol26/projects/wssd_wgac/nobackups/dgordon/wgac/panTro6/fasta/chrUn_NW_019932986v1:0:451 -j /net/eichler/vol26/projects/wssd_wgac/nobackups/dgordon/wgac/panTro6/fasta/chrUn_NW_019933313v1:1231:1231 -l 110000 -o /net/eichler/vol26/projects/wssd_wgac/nobackups/dgordon/wgac/panTro6/data/align_both/0022/both0113351 -f -40 -g -1 -b 400000 && cd /net/eichler/vol26/projects/wssd_wgac/nobackups/dgordon/wgac/panTro6 && rm -rf /net/eichler/vol26/projects/wssd_wgac/nobackups/dgordon/wgac/panTro6/both_tmp/both0113351 && touch global_align_110k_flags/11
                    szCommand = szLine
                    aCommandWords = szCommand.split()
                    bFoundMkdir = True

                    fSummarizeOutput.write( "command: " + szCommand )
                    fSummarizeOutput.write( "error lines: " )

                elif( bFoundMkdir ):
                    nLinesAfterCommand += 1

                    szErrorLines += " "
                    szErrorLines += szLine.rstrip()
                    
                    if ( nLinesAfterCommand >= 10 ):
                        break
            # for szLine in fLog.readlines():

            fSummarizeOutput.write( szErrorLines + "\n" )
            if ( "BAD MPOS" in szErrorLines ):
                nBadMpos += 1
            else:
                nNoBadMpos += 1
            
    #for szJobID in aJobIDs:
    
    fSummarizeOutput.write( "BAD MPOS: {:d} not BAD MPOS: {:d}\n".format( nBadMpos, nNoBadMpos ) )

    szCommand = "wc -l sge_align_batch.sh"
    print "about to execute: " + szCommand
    szOutput = subprocess.check_output( szCommand, shell = True )

    nSgeAlignBatchLines = int( szOutput.split()[0] )

    fPerCentBad = float( nBadMpos + nNoBadMpos ) * 100.0 / nSgeAlignBatchLines

    fSummarizeOutput.write( szOutput )
    fSummarizeOutput.write( "error jobs = {:.2f} %\n".format( fPerCentBad ) )


    
print "\n\n\nSee: summarize_global_align_errors.txt\n\n"


    
    



