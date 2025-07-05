#!/bin/bash -l
#SBATCH --time=36:00:00
#SBATCH --ntasks=128
#SBATCH --mem=200gb
#SBATCH --nodes=1
#SBATCH --mail-type=ALL
#SBATCH --mail-user=shoops@umn.edu
#SBATCH --partition=amd512
#SBATCH --account=knightsd

source ~/.bashrc

cd /scratch.global/shoops/

# Running metaSPAdes on Donor and Pre samples
## donors (closed label)
for donor in PREP_2019_00017 PREP_2019_00018 PREP_2019_00019 PREP_2021_00012 PREP_2021_00117
do
echo "Run metaSPADes for donor ${donor}..."
mkdir metaspades/${donor}
~/lib/SPAdes-3.15.5-Linux/bin/metaspades.py -t 128 -1 ucfmt_kneaddata/${donor}.r1_kneaddata_paired_1.fastq -2 ucfmt_kneaddata/${donor}.r1_kneaddata_paired_2.fastq -o metaspades/${donor}/
done
## pre samples (closed label)
for pre in UCFMT001_BL UCFMT004_BL UCFMT017_BL UCFMT023_BL UCFMT026_BL UCFMT028_BL
do
echo "Run metaSPADes for pre ${pre}..."
mkdir metaspades/${pre}
~/lib/SPAdes-3.15.5-Linux/bin/metaspades.py -t 128 -1 ucfmt_kneaddata/${pre}.r1_kneaddata_paired_1.fastq -2 ucfmt_kneaddata/${pre}.r1_kneaddata_paired_2.fastq -o metaspades/${pre}/
done

echo "Done."

# Bowtie2-build commands
bowtie2-build --seed 25 -threads 128 metaspades/${donor}/scaffolds.fasta bowtie_db/${donor}/


# Bowtie2-align commands
time bowtie2-align-s -x ./bowtie2_db/${donor} -S ./bowtie2_out/${post}_donor_alignment_1.00.txt --very-sensitive --np 1 --mp "1,1" --rdg "0,1" --rfg "0,1" --score-min "C,0" -k 1  -p 16 --no-hd --no-unal -q ${post}.r1.kneaddata_paired_1.fastq


