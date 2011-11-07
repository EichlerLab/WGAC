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

scriptName='/net/eichler/vol2/local/inhousebin/align_fast2_batch4.pl'
scriptParams='-u both:random -n both -d 5000 -i 0:1:2 -j 4:5:6 -t both.parse.defugu.trim.fixed.trim.defrac -o align_both';


# translate 1st arg to full path
fastaDir=$(readlink -f $1)

# translate 2nd arg
if [ "$2" == "64bit" ]
then
	clusterName='eeek'
elif [ "$2" == "32bit" ]
then
	clusterName='eee'
else
	clusterName=''
fi

if [ -z "$clusterName" ]
then
	# 
	perl $scriptName $scriptParams -l 110000 -f $fastaDir 
else
	# cd to where the both_tmp directory will be created, then run Perl 
	echo "using nodes of $clusterName to run $scriptName ..."
	cwd=$(pwd)
	cexec $clusterName: "cd $cwd; perl $scriptName $scriptParams -f $fastaDir >/dev/null" >/dev/null
fi

 



