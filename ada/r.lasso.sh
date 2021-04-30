#BSUB -L /bin/bash              # uses the bash login shell to initialize the job's execution environment.
#BSUB -J R.lasso           # job name
#BSUB -n 1                     # assigns 10 cores for execution
#BSUB -R "span[ptile=1]"       # assigns 10 cores per node
#BSUB -R "rusage[mem=256000]"     # reserves 1000MB memory per core
#BSUB -M 256000                   # sets to 1,000MB (~1GB) per process enforceable memory limit. Total job memory = (M * n)
#BSUB -W 24:00                   # sets to 96 hours the job's runtime wall-clock limit.
#BSUB -o /scratch/user/simonpan/ct/routput/R.lasso.log     # directs the job's standard output to stdout.jobid
#BSUB -e /scratch/user/simonpan/ct/routput/R.lasso.error.log     # directs the job's standard error to stderr.jobid
#BSUB -u simonpan@exchange.tamu.edu       #Send all emails to email_address
#BSUB -B -N

module load R/4.0.3-foss-2020a-recommended-mt

Rscript /scratch/user/simonpan/ct/quafunreg_ada.R

#  cd /scratch/user/simonpan/MOSJTumor/temp/tmp && bsub < /scratch/user/simonpan/ct/r.lasso.sh