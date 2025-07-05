#!/bin/bash -l
#SBATCH --time=48:00:00
#SBATCH --ntasks=128
#SBATCH --mem=200gb
#SBATCH --nodes=1
#SBATCH --mail-type=ALL
#SBATCH --mail-user=shoops@umn.edu
#SBATCH --partition=amd512
#SBATCH --account=knightsd

source ~/.bashrc
cd /scratch.global/shoops/
module load bowtie2


##### METASPADES - PLACEBO PRE SAMPLES #####
echo "Starting metaSPAdes for donor groups..."
## NEEDS FILE: placebo_pre_table.txt
preinput="./placebo_pre_table.txt"
i=1
while ICF=$'\t' read -r -a myArray
do
  # skip first line (header)
  test $i -eq 1 && ((i=i+1)) && continue
  # run metaspades on respective forward & reverse reads (concatenated files)
  echo "Run metaSPADes for pre ${myArray[0]}..."
  mkdir metaspades/pre_${myArray[0]}
  ~/lib/SPAdes-3.15.5-Linux/bin/metaspades.py -t 128 -1 kneaddata_out/${myArray[1]}_kneaddata_paired_1.fastq -2 kneaddata_out/${myArray[1]}_kneaddata_paired_2.fastq -o metaspades/pre_${myArray[0]}
done < "$preinput"
echo "Completed metaSPAdes for placebo pre samples...\n"


##### BOWTIE2-BUILD PLACEBO PRE SAMPLES #####
echo "Starting bowtie2-build for donor groups..."
i=1
while ICF=$'\t' read -r -a myArray
do
  # skip first line (header)
  test $i -eq 1 && ((i=i+1)) && continue
  # build bowtie2 databases for placebo pre samples
  bowtie2-build --seed 25 metaspades/pre_${myArray[0]}/scaffolds.fasta bowtie_db/pre_${myArray[0]}
done < "$preinput"
chmod -R 750 bowtie_db/
echo "Completed bowtie2-build for placebo pre samples.\n"


##### BOWTIE2-ALIGN PLACEBO POST SAMPLES #####
## NEED FILE : post_table.txt
postinput="./placebo_post_table.txt"
## note: post table is organized as : ID | Post file ID | Donor | Donor_Set | Pre (0-indexing)
i=1
while IFS=$'\t' read -r -a myArray
do
  # skip first line (header)
  test $i -eq 1 && ((i=i+1)) && continue
  # bowtie2 alignment per post sample - 100% match
  /usr/bin/time bowtie2-align-s -x ./bowtie_db/donor_${myArray[3]} -S ./bowtie_out_placebo/${myArray[1]}_donor_${myArray[3]}_alignment_1.00.txt --very-sensitive --np 1 --mp "1,1" --rdg "0,1" --rfg "0,1" --score-min "C,0" -k 1  -p 16 --no-hd --no-unal -q ./kneaddata_out/${myArray[1]}_kneaddata_paired_1.fastq
  /usr/bin/time bowtie2-align-s -x ./bowtie_db/pre_${myArray[0]} -S ./bowtie_out_placebo/${myArray[1]}_pre_${myArray[0]}_alignment_1.00.txt --very-sensitive --np 1 --mp "1,1" --rdg "0,1" --rfg "0,1" --score-min "C,0" -k 1  -p 16 --no-hd --no-unal -q ./kneaddata_out/${myArray[1]}_kneaddata_paired_1.fastq
  # bowtie2 alignment per post sample - 99% match
  /usr/bin/time bowtie2-align-s -x ./bowtie_db/donor_${myArray[3]} -S ./bowtie_out_placebo/${myArray[1]}_donor_${myArray[3]}_alignment_0.99.txt --very-sensitive --np 1 --mp "1,1" --rdg "0,1" --rfg "0,1" --score-min "L,0,-.01" -k 1  -p 16 --no-hd --no-unal -q ./kneaddata_out/${myArray[1]}_kneaddata_paired_1.fastq
  /usr/bin/time bowtie2-align-s -x ./bowtie_db/pre_${myArray[0]} -S ./bowtie_out_placebo/${myArray[1]}_pre_${myArray[0]}_alignment_0.99.txt --very-sensitive --np 1 --mp "1,1" --rdg "0,1" --rfg "0,1" --score-min "L,0,-.01" -k 1  -p 16 --no-hd --no-unal -q ./kneaddata_out/${myArray[1]}_kneaddata_paired_1.fastq
  # bowtie2 alignment per post sample - 98% match
  /usr/bin/time bowtie2-align-s -x ./bowtie_db/donor_${myArray[3]} -S ./bowtie_out_placebo/${myArray[1]}_donor_${myArray[3]}_alignment_0.98.txt --very-sensitive --np 1 --mp "1,1" --rdg "0,1" --rfg "0,1" --score-min "L,0,-.02" -k 1  -p 16 --no-hd --no-unal -q ./kneaddata_out/${myArray[1]}_kneaddata_paired_1.fastq
  /usr/bin/time bowtie2-align-s -x ./bowtie_db/pre_${myArray[0]} -S ./bowtie_out_placebo/${myArray[1]}_pre_${myArray[0]}_alignment_0.98.txt --very-sensitive --np 1 --mp "1,1" --rdg "0,1" --rfg "0,1" --score-min "L,0,-.02" -k 1  -p 16 --no-hd --no-unal -q ./kneaddata_out/${myArray[1]}_kneaddata_paired_1.fastq
done < "$postinput"
echo "Completed bowtie2-align with placebo samples.\n"


##### DONE #####
echo 'ALL DONE.'