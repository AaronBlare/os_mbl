# !/bin/bash -i
#$ -S /bin/bash

# --- Mandatory qsub arguments
# Hardware requirements.
#$ -l h_rss=4000M,h_fsize=4000M,h_cpu=80:00:00,hw=x86_64

# --- Optional qsub arguments
# Change working directory - your job will be run from the directory
# that you call qsub in. So stdout and stderr will end up there .
# $ -cwd

data_path=/data/biophys/yusipov/mbl_zero_super_vector
scratch=/scratch/yusipov/mbl_zero_super_vector/$1
code_base=$HOME/mbl_zero_super_vector/source/
mkdir -p $scratch
mkdir -p $1
cd $scratch
cp $code_base/* .
cp $1/config.txt .

echo " Running on $(hostname)"
echo " We are in $(pwd) "

cat config.txt
matlab -nodesktop -nosplash -r mbl_zero_super_vector

cp -r $scratch/*.txt $1 
rm -r $scratch/*


