#! /bin/sh

# Usage: qsub -q all.q repeatOut.sh

# Specify the shell for this job
#$ -S /bin/bash 
export MPICH_PROCESS_GROUP=no
export P4_RSHCOMMAND=/usr/bin/rsh
export PRINT_SEQUENCES=2
#export C3_RSH='ssh -q'
#ulimit -c 0

# pe request
#$ -pe mpich 10-15

# Send an email when the script begins, ends, aborts, or suspends.
#$ -m beas

#$ -hard 

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

# Set defaults.
WORKING_DIR=$SGE_O_WORKDIR
INPUT=$WORKING_DIR/fasta
OUTPUT=$WORKING_DIR/mask_out

# Get options from the user.
while getopts :i:o: OPTION
do
  case $OPTION in
    i)
      INPUT=$OPTARG
      ;;
    o)
      OUTPUT=$OPTARG
      ;;
    \?)
      echo "Invalid option: -$OPTARG" >&2
      exit 1
      ;;
    :)
      echo "Option -$OPTARG requires an argument." >&2
      exit 1
      ;;      
  esac
done

/net/eichler/vol6/software/mpich-1.2.7-amd64/bin/mpirun -np $NSLOTS -machinefile $TMPDIR/machines /net/eichler/vol5/home/ssajjadi/wgacbin/step_2_mpi/repeatout "perl /net/eichler/vol5/home/tinlouie/wgacscripts/maskOutGenFromLowCase.pl"  $INPUT $OUTPUT $JOB_ID


