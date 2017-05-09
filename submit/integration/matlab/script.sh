# !/bin/bash -i
#$ -S /bin/bash

# --- Mandatory qsub arguments
# Hardware requirements.
#$ -l h_rss=4000M,h_fsize=4000M,h_cpu=80:00:00,hw=x86_64

# --- Optional qsub arguments
# Change working directory - your job will be run from the directory
# that you call qsub in. So stdout and stderr will end up there .
# $ -cwd

# --- Job Execution
# For faster disk access copy files to / scratch first .
module load mkl

scratch=/scratch/yusipov/matlab/$1
my_folder=$HOME/MBL/matlab/integration
code_base=$HOME/MBL/matlab/integration/sources
mkdir -p $scratch
mkdir -p $my_folder/$1
cd $scratch
cp $code_base/* .
cp $my_folder/$1/config.txt .

# Execution - running the actual program .
# [ Remember : Don ’ t read or write to / home from here .]

echo " Running on $(hostname)"
echo " We are in $(pwd) "

cat config.txt
matlab -nodesktop -nosplash -r mbl_single

# Finish - Copy files back to your home directory , clean up .
cp -r $scratch/*.txt $my_folder/$1 # Better use a subdirectory of $HOME .
rm -r $scratch/*
#cd
#rm - rf $scratch
