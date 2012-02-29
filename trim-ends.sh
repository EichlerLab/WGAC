#!/bin/bash

# Trim ends of WGAC sequences.
#
# Usage: qsub -q all.q trim-ends.sh -i <INPUT> -f <FASTA FILES> -o <OUTPUT>

# Specify the shell for this job
#$ -S /bin/bash 

export MPICH_PROCESS_GROUP=no
export P4_RSHCOMMAND=/usr/bin/rsh
export PRINT_SEQUENCES=2
export C3_RSH='rsh -q'
ulimit -c 0

export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/usr/lib64/openmpi/1.4-gcc/lib
export PATH=$PATH:/usr/lib64/openmpi/1.4-gcc/bin
export MANPATH=$MANPATH:/usr/lib64/openmpi/1.4-gcc/man

# pe request
#$ -pe orte 20-50

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

#First agrument is the command
#command="\"perl /net/eichler/vol5/home/ssajjadi/wgacbin/step_8_mpi/Trim.pl\"";

# Set defaults.
WORKING_DIR=$SGE_O_WORKDIR
INPUT=$WORKING_DIR/data/step_8_mpi/defugu
FASTA=$WORKING_DIR/fasta
OUTPUT=$WORKING_DIR/data/step_8_mpi/trim

# Get options from the user.
while getopts :i:f:o: OPTION
do
  case $OPTION in
    i)
      INPUT=$OPTARG
      ;;
    f)
      FASTA=$OPTARG
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

mpirun -mca btl ^openib -np $NSLOTS \
  /net/eichler/vol4/home/jlhudd/wgac/trim_ends/runTrim \
  "perl /net/eichler/vol4/home/jlhudd/wgac/trim_ends/Trim.pl" \
  $INPUT $FASTA $OUTPUT $JOB_ID

