#!/bin/bash
# 1)  THIS SCRIPT WILL RUN cexec, SO YOU MUST RUN THIS SCRIPT ON eeek:
# 2)  THIS SCRIPT ASSUMES THAT PREVIOUS STEP OF WGAC PIPELINE PRODUCED A FILE CALLED both.parse.defugu.trim.fixed.trim.defrac
#     & IT IS IN THE CURRENT WORKING DIRECTORY

# usage:
# wgacAlign.sh ../fasta 64bit
# wgacAlign.sh ../fasta
#
# 1st arg:  directory containing fractionated fasta files (e.g. /net/eichler/vol2/eee_shared/GRCh37/fasta)
# 2nd arg:  if the arg is the phrase 64bit, then cexec will use the 64-bit cluster (eeek)
#           if                       32bit,                         32-bit cluster (eee)
#           otherwise, skip cexec and run Perl program on single node to wrap up loose ends

scriptName='/net/eichler/vol2/home/dgordon/wgac/github/WGAC/align_fast2_batch4.pl'

scriptParams='-u both:random -n both -d 5000 -b 400000 -i 0:1:2 -j 4:5:6 -t data/both.parse.defugu.trim.fixed.trim.defrac -o data/align_both';

# translate 1st arg to full patih

fastaDir=$(readlink -f $1)

align_len=$2
	
perl $scriptName $scriptParams -f $fastaDir -l $align_len

