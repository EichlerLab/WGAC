#! /bin/bash
#$ -S /bin/bash
#$ -cwd
#$ -V

export P4_RSHCOMMAND=/usr/bin/rsh
export PRINT_SEQUENCES=2
export JOBDIR=$TMPDIR/work
export WORKING_DIR=$SGE_O_WORKDIR
export PATH=$PATH:/net/eichler/vol2/local/bin

if [[ "$#" -lt "1" ]]
then
    echo "Usage: $0 <command list> [env.sh]"
    exit 1
fi

INPUT=$1
ENV_FILE=$2

. /etc/profile.d/modules.sh

if test ! -z $MODULESHOME; then
   module load modules modules-init/prod modules-gs/prod
   module load modules-eichler/prod
   module load openmpi/1.5.4
fi

module unload python
module load python/2.7.2

if [[ ! -z ${ENV_FILE} ]]
then
    source ${ENV_FILE}
    echo "Loaded environment from ${ENV_FILE}"
fi

if [[ ! -e ${INPUT} ]]
then
    echo "Command list file doesn't exist: ${INPUT}"
    exit 1
fi

echo "Python:"
which python
echo "Python path:"
echo $PYTHONPATH
echo `hostname`

mpirun -x PATH -x LD_LIBRARY_PATH \
  --prefix $MPIBASE -mca plm ^rshd \
  -mca btl ^openib \
  /net/eichler/vol4/home/jlhudd/src/general_pipe/general_pipe.py \
  --input_file="$INPUT" \
  --commands
echo `hostname`
