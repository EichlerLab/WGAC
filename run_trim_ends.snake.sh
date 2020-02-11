module purge
module load modules modules-init modules-gs/prod modules-eichler/prod libevent/2.1.3-alpha

module load miniconda/4.5.12



snakemake -s trim_ends.snake  --jobname "{rulename}.{jobid}" --drmaa " -V -cwd -e ./log -o ./log {params.sge_opts}  -S /bin/bash" -w 100 --jobs 100 -p -k
