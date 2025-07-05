##### RENAMING & KNEADDATA PREP #####

# get the database :
#kneaddata_database --download human_genome bowtie2 human_genome_db
# paired end kneaddata format :
#kneaddata --input1 $INPUT1 --input2 $INPUT2 --reference-db $DATABASE --output $OUTPUT_DIR

# rename the identifiers so they match for kneaddata
cd /scratch.global/shoops/raw_ucfmt/
for fn1 in *R1_001.fastq
do
  sed 's/ 1:N:0.*$/\/1/g' < $fn1 > ../raw_ucfmt_rename/$fn1
done
for fn2 in *R2_001.fastq
do
  sed 's/ 2:N:0.*$/\/2/g' < $fn2 > ../raw_ucfmt_rename/$fn2
done



##### SLURM SCRIPT BELOW #####

#!/bin/bash -l
#SBATCH --time=72:00:00
#SBATCH --ntasks=128
#SBATCH --mem=120gb
#SBATCH --nodes=1
#SBATCH --mail-type=ALL
#SBATCH --mail-user=shoops@umn.edu
#SBATCH --partition=amdsmall
#SBATCH --account=knightsd

cd /scratch.global/shoops/

# load modules
source ~/.bashrc
conda activate biobakery_workflow

# conduct kneaddata for a single file to test - WORKED! Took about 3 hrs...
#echo "TESTING: Kneaddata for one file"
#time kneaddata --input1 raw_ucfmt/UCFMT001BL_S55_R1_001.fastq.gz --input2 raw_ucfmt/UCFMT001BL_S55_R2_001.fastq.gz -t 128 --reference-db human_genome_db --output test_out

# conduct kneaddata for all read files
echo "START OF LOOP"
for fn in raw_ucfmt_rename/*_R1_001.fastq.gz
do
  time kneaddata -t 128 --input1 $fn --input2 ${fn//_R1_/_R2_} --reference-db human_genome_db --output kneaddata_out
done

# close conda environment
conda deactivate
echo "ALL DONE.


