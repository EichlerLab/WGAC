#! /bin/sh
# Specify the shell for this job
#$ -S /bin/sh 
export MPICH_PROCESS_GROUP=no
export P4_RSHCOMMAND=/usr/bin/rsh
export PRINT_SEQUENCES=2
export C3_RSH='rsh -q'
ulimit -c 0
 # Tell Sun Grid Engine to send an email when the job begins
 # and when it ends.
# pe request
#$ -pe mpich 3
#$ -M ssajjadi@u.washington.edu
#send it when it started&finished
#$ -m beas
#$ -hard 
# Specify the location of the output
#$ -o /net/eichler/vol5/home/ssajjadi/step_8_13/test
#$ -e /net/eichler/vol5/home/ssajjadi/step_8_13/test
# Location of executables
progpath=/mnt/local/bin
echo "Got $NSLOTS slots" 
echo "path=$PATH" 
echo "P4_RSHCOMMAND=$P4_RSHCOMMAND" 
echo "machine_file=$TMPDIR/machines" 
echo "JOB_ID=$JOB_ID" 
echo "TEMDPIR=$TMPDIR" 
echo "HOSTNAME=$HOSTNAME"
#sort -u -o $TMPDIR/Mfile $TMPDIR/machines
#echo $HOSTNAME >> $TMPDIR/Mfile
#thisslots=`wc -l $TMPDIR/Mfile|awk '{print $1}'`
#echo "thisslots=$thisslots"
#4std; mean: 209.57; std: 13.40

#First agrument is the command
#command="\"perl /net/eichler/vol5/home/ssajjadi/wgacbin/step_8_mpi/Trim.pl\"";

#Second argument is the location of input files (defugu files)
inputdir="/net/eichler/vol5/home/ssajjadi/step_8_13/test/defugu";

#Third argument is the location of fasta files
fastadir="/net/eichler/vol5/home/ssajjadi/fasta";

#Fourth argument is the location of output files (trim files)
outputdir="/net/eichler/vol5/home/ssajjadi/step_8_13/test/trim";

/net/eichler/vol6/software/mpich-1.2.7-amd64/bin/mpirun -np $NSLOTS -machinefile $TMPDIR/machines /net/eichler/vol5/home/ssajjadi/wgacbin/step_8_mpi/runTrim "perl /net/eichler/vol5/home/ssajjadi/wgacbin/step_8_mpi/Trim.pl"  $inputdir $fastadir $outputdir $JOB_ID

