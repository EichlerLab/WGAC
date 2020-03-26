mkdir -p log
module purge
module load modules modules-init modules-gs/prod modules-eichler/prod miniconda/4.5.12
snakemake -s lastz_self.snake --jobname "{rulename}.{jobid}" --drmaa " -V -cwd -e ./log -o ./log {params.sge_opts} -l testing=true  -S /bin/bash" -w 100 --jobs 300 -p -k

