# !/bin/bash -i
#$ -S /bin/bash

# --- Mandatory qsub arguments
# Hardware requirements.
#$ -l h_rss=4000M,h_fsize=40M,h_cpu=10:00:00,hw=x86_64

# --- Optional qsub arguments
# Change working directory - your job will be run from the directory
# that you call qsub in. So stdout and stderr will end up there .
# $ -cwd

# --- Job Execution
# For faster disk access copy files to / scratch first .

scratch=/scratch/yusipov
rm -r $scratch/*
