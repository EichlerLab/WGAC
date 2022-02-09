module purge

module load modules modules-init modules-gs/prod modules-eichler/prod

module load miniconda/4.5.12
snakemake -s trf.snake populate_trf_with_links

# added Feb 9, 2022 (DG):
# restart-times was necessary since occasionally got this error:
# Error recording metadata for finished job ([Errno 2] No such file or directory: 'fugu_trf/chr2_274.fugu'). Please ensure write permissions for the directory /net/eichler/vol27/projects/hprc/nobackups/chm13v2_wgac3/.snakemake

snakemake -s trf.snake --jobname "{rulename}.{jobid}"  --drmaa " -w n -V -cwd -w n -e ./log -o ./log {params.sge_opts} -S /bin/bash" -j 100 -k --rerun --restart-times 1
