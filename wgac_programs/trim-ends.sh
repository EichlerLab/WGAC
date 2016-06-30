#!/bin/bash

# Trim ends of WGAC sequences.
#
# Usage: qsub -q all.q trim-ends.sh -i <INPUT> -f <FASTA FILES> -o <OUTPUT>

# Specify the shell for this job
#$ -S /bin/bash

set -e

. /etc/profile.d/modules.sh

if test ! -z $MODULESHOME; then
   module load modules modules-init modules-gs modules-eichler
   module load openmpi/1.5.3
fi

module load python/2.7.2

export P4_RSHCOMMAND=/usr/bin/rsh
export PRINT_SEQUENCES=2
export C3_RSH='rsh -q'
ulimit -c 0

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

# Create tmp directory on cluster nodes and make sure it's empty.
TRIM_TMP_DIR=/tmp/endtrim

mpirun -x PATH -x LD_LIBRARY_PATH \
  --prefix $MPIBASE -mca plm ^rshd \
  -mca btl ^openib python /net/eichler/vol7/home/psudmant/EEE_Lab/projects/batch_node_copy/code/batch_node_copy.py \
  --source "~jlhudd/wgac/trim-ends.sh" --dest "/dev/null" \
  --pre_sync_commands "rm -rf ${TRIM_TMP_DIR}; mkdir -p ${TRIM_TMP_DIR}"

echo "Trimming ends"
mpirun -x PATH -x LD_LIBRARY_PATH \
  --prefix $MPIBASE -mca plm ^rshd \
  -mca btl ^openib \
  /net/eichler/vol4/home/jlhudd/wgac/trim_ends/runTrim \
  "perl /net/eichler/vol4/home/jlhudd/wgac/trim_ends/Trim.pl" \
  $INPUT $FASTA $OUTPUT $JOB_ID

