#! /bin/sh

# Run BLAST64 for the given input sequences.
#
# Usage: qsub -q all.q blast64.sh -i <INPUT> -o <OUTPUT> -b <BLASTDB>

# Specify the shell for this job
#$ -S /bin/sh

export MPICH_PROCESS_GROUP=no 
export P4_RSHCOMMAND=/usr/bin/rsh 
export JOBDIR=$TMPDIR/work

#$ -pe mpich 40

# Send an email when the script begins, ends, aborts, or suspends.
#$ -m beas

# The job is located in the current
# working directory.
#$ -cwd

echo "Got $NSLOTS slots"
echo "path=$PATH"
echo "P4_RSHCOMMAND=$P4_RSHCOMMAND"
echo "machine_file=$TMPDIR/machines"

WORKING_DIR=$SGE_O_WORKDIR
INPUT=$WORKING_DIR/fugu
OUTPUT=$WORKING_DIR/blastout
DATABASE_PATH=$WORKING_DIR/blastdb/bofugu
BLAST_TMP_DIR=/var/tmp/blastdb

# Get options from the user.
while getopts :i:o:b:t: OPTION
do
  case $OPTION in
    i)
      INPUT=$OPTARG
      ;;
    o)
      OUTPUT=$OPTARG
      ;;
    b)
      DATABASE_PATH=$OPTARG
      ;;
    t)
      BLAST_TMP_DIR=$OPTARG
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

# Get short name of the database.
DATABASE_NAME=`basename $DATABASE_PATH`

echo Database name: $DATABASE_NAME
echo Database path: $DATABASE_PATH
echo Tmp path: $BLAST_TMP_DIR
echo Input: $INPUT
echo Output: $OUTPUT

# Create tmp directory on cluster nodes.
/usr/bin/uniq $TMPDIR/machines | /net/eichler/vol2/local/bin/rgang - mkdir -p $BLAST_TMP_DIR

# Copy BlastDB files to nodes.
/usr/bin/uniq $TMPDIR/machines | /net/eichler/vol2/local/bin/rgang --rcpt 8000  -c --rcp="/usr/bin/rcp -r"  - $DATABASE_PATH.[^q]* $BLAST_TMP_DIR

# Run BLAST.
#/net/eichler/vol6/software/mpich-1.2.7-i386/bin/mpirun -np $NSLOTS -machinefile $TMPDIR/machines /net/eichler/vol7/home/ginger/bin/general_pipe "/net/eichler/vol2/local/bin/blastall -p blastn -i dummy_in -o dummy_out  -d $BLAST_TMP_DIR/$DATABASE_NAME -G 180 -E 1 -q -80 -r 30 -e 1e-30 -W 28"   $INPUT    $OUTPUT ".bo" $JOBDIR

/net/eichler/vol6/software/mpich-1.2.7-i386/bin/mpirun -np $NSLOTS -machinefile $TMPDIR/machines /net/eichler/vol7/home/ginger/bin/general_pipe "/net/eichler/vol2/local/bin/blastall -p blastn -i dummy_in -o dummy_out  -d $BLAST_TMP_DIR/$DATABASE_NAME -z 3000000000 -Y 3000000000 -v 5000 -b 5000 -G 180 -E 1 -q -80 -r 30 -e 1e-30 -F F"   $INPUT    $OUTPUT ".bo" $JOBDIR

# Create tmp directory on cluster nodes.
/usr/bin/uniq $TMPDIR/machines | /net/eichler/vol2/local/bin/rgang - rm -rf $BLAST_TMP_DIR


# blastall is v2.2.11, not much different than 2.2.20
# -p blastn (program name)
# -d /net/eichler/vol4/home/linchen2/wgac/GRCh37/blastdb/bofugu (database)

# -z 3000000000 (?? effective length of database, useful for maintaining consistent statistics as databases grow ??)
# -Y 3000000000 (?? effective length of the search space = db size * query ??)
# http://etutorials.org/Misc/blast/Part+V+BLAST+Reference/Chapter+13.+NCBI-BLAST+Reference/13.3+blastall+Parameters/

# -v 5000 (?? # of database seqs to show one-line descriptions)
# -b 5000 (# of database seqs to show alignments)

# -G 180 (gap open cost)
# -E 1 (gap extension cost)
# -q -80 (nuc mismatch penalty)
# -r 30 (nuc match reward)
# -e 1e-30 (e-value)
# -F F (don't filter query seq)
