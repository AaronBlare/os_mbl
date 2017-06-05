# !/bin/bash -i
#$ -S /bin/bash

# --- Mandatory qsub arguments
# Hardware requirements.
#$ -l h_rss=4000M,h_fsize=4000M,h_cpu=48:00:00,hw=x86_64

# --- Optional qsub arguments
# Change working directory - your job will be run from the directory
# that you call qsub in. So stdout and stderr will end up there .
# $ -cwd

# --- Job Execution
# For faster disk access copy files to / scratch first .
module load mkl

scratch=/scratch/yusipov/matlab/os_mbl/$1
my_folder=/data/biophys/yusipov/os_mbl/int/matlab
code_base=$HOME/Work/os_mbl/matlab/plot/integration/imbalance
mkdir -p $scratch
mkdir -p $1
cd $scratch
cp $code_base/imbalance_pdf_dump.m .
cp $1/config.txt .

# Execution - running the actual program .
# [ Remember : Don ’ t read or write to / home from here .]

echo " Running on $(hostname)"
echo " We are in $(pwd) "

cat config.txt

matlab -nodesktop -nosplash -r imbalance_pdf_dump

# Finish - Copy files back to your home directory , clean up .
cp -r $scratch/*.txt $1 # Better use a subdirectory of $HOME .
rm -r $scratch/*
#cd
#rm - rf $scratch
