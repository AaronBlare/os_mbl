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
module load mkl

scratch=/scratch/yusipov/MBL_TR/$1
code_base=$HOME/mbl_transition_rates/cpp/source
mkdir -p $scratch
mkdir -p $HOME/mbl_transition_rates/cpp/$1
cd $scratch
cp $2/$1/config.txt .

# Execution - running the actual program .
# [ Remember : Don ’ t read or write to / home from here .]

echo " Running on $(hostname)"
echo " We are in $(pwd) "

cat config.txt
$code_base/tr.out

# Finish - Copy files back to your home directory , clean up .
cp -r $scratch/* $2/$1 # Better use a subdirectory of $HOME .
rm -r $scratch/*
#cd
#rm - rf $scratch
