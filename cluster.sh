# (1) Run the job in the same directory that qsub is issued
#
#$ -cwd
#
# (2) Give the job a meaningful name (this will appear in qstat listings)
#
#$ -N sound_speed_analysis
#
# (3) Send a notification email when your job begins (b) and ends (e)
# being run
#
#$ -m be
#
# (4) Set the queue that the job should run on. As undergraduates, you have
# access to the sopa.1.day queue.
#
#$ -q sopa.1.day
#
# (5) Give an estimate of how long the job should run for, expressed as
# hours:minutes:seconds. This can be up to 24 hours on the sopa.1.day queue.
# NB: your job will be halted when this time expires, so it is best to slightly
# over-estimate.
#
#$ -l h_rt=04:00:00

#$ -pe omp 16

export OMP_NUM_THREADS=${NSLOTS:-4}
build/ProjectModelling ${1:-data}
 
# The ${n:-default} syntax looks for the n'th argument supplied to the
# script, and substitutes default if this is empty.
