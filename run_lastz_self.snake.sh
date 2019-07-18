mkdir -p log
module purge
module load modules modules-init modules-gs/prod modules-eichler/prod miniconda/4.5.12
snakemake -s lastz_self.snake --drmaa " -q eichler-short.q -l h_rt=24:00:00  -V -cwd -e ./log -o ./log {params.sge_opts}  -S /bin/bash" -w 100 --jobs 100 -p -k
