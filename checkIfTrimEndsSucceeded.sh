nstep_8_mpi_defugu=`ls data/step_8_mpi/defugu/* | wc -l`
nstep_8_mpi_trim_size=`ls data/step_8_mpi/trim/*/* | wc -l`

if [ $nstep_8_mpi_defugu -eq $nstep_8_mpi_trim_size ]
then
    echo "data/step_8_mpi/defugu and data/step_8_mpi/trim are the same size"
    exit 0
else
    echo "data/step_8_mpi/defugu size  $nstep_8_mpi_defugu"
    echo "data/step_8_mpi/trim size $nstep_8_mpi_trim_size"
    echo "data/step_8_mpi/defugu and data/step_8_mpi/trim are not the same size so qsub -N trim-ends must have failed to complete successfully"
    exit 1
fi
