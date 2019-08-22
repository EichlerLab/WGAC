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

# changed (DG, 8/16/2019) so Trim.pl is part of wgac distribution instead of in jlhudd's directory
for i in `pwd`/data/step_8_mpi/defugu/*;do echo "perl `pwd`/Trim.pl $i `pwd`/fasta `pwd`/data/step_8_mpi/trim ${SPECIES}";done

