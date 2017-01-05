#!/bin/bash

set -e

if [[ "$#" -ne "1" ]]
then
    echo "Usage: $0 <species>"
    exit 1
fi

SPECIES=$1
echo "#!/bin/bash"
echo "#$ -S /bin/bash"
echo "#$ -cwd"

for i in `pwd`/data/step_8_mpi/defugu/*;do echo "perl /net/eichler/vol4/home/jlhudd/wgac/Trim.pl $i `pwd`/fasta `pwd`/data/step_8_mpi/trim ${SPECIES}";done

