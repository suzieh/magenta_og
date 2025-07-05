#!/bin/bash -l
#SBATCH --time=84:00:00
#SBATCH --ntasks=128
#SBATCH --mem=250gb
#SBATCH --nodes=1
#SBATCH --mail-type=ALL
#SBATCH --mail-user=shoops@umn.edu
#SBATCH --partition=amd512
#SBATCH --account=knightsd

# Script for running additional simulations
cd /scratch.global/shoops/smillie
source /home/knightsd/shoops/.bashrc

# --- Clean up reads ---
echo "--- Begin by QC cleaning/trimming the FASTQs...."
cd PRJEB23524/
conda activate biobakery_workflow
for fn in *_1.fastq.gz
do
  time kneaddata -t 128 --input1 $fn --input2 ${fn//_1/_2} --reference-db human_genome_db --output kneaded
done
conda deactivate
# note kneaddata output we want: kneaded/ERR#######_kneaddata_paired_1.fastq
echo "....Completed QC."


# --- Make mixtures per donor/pre pair ---
# amounts : 0 875000 1750000 2625000 3500000 4375000 5250000 6125000 7000000 7875000 8750000
echo "--- Generating mixtures for output files...."
cd /scratch.global/shoops/smillie
tinput="./triads_table.txt"
i=1
while ICF=$'\t' read -r -a myArray
do
  # skip header for now
  test $i -eq 1 && ((i=i+1)) && continue
  # get donor ascension and pre ascension
  echo "Mixing ${myArray[0]}: Donor ${myArray[9]} & Pre ${myArray[8]}"
  # make 0% donor mixture (forward and reverse)
  head -n 8750000 PRJEB23524/kneaded/${myArray[8]}_kneaddata_paired_1.fastq > ./mixtures/${myArray[0]}_simulation_0.fastq
  # make 100% donor mixture
  head -n 8750000 PRJEB23524/kneaded/${myArray[9]}_kneaddata_paired_1.fastq > ./mixtures/${myArray[0]}_simulation_100.fastq
  # loop over remaining mixtures (forward only)
  for p in 10 20 30 40 50 60 70 80 90
  do
    # add donor mixture
    donor_share=`expr 87500 \\* $p`
    head -n $donor_share PRJEB23524/kneaded/${myArray[9]}_kneaddata_paired_1.fastq > ./mixtures/${myArray[0]}_simulation_$p.fastq
    # add pre mixture
    invp=`expr 100 - $p`
    pre_share=`expr 87500 \\* $invp`
    head -n $pre_share PRJEB23524/kneaded/${myArray[8]}_kneaddata_paired_1.fastq >> ./mixtures/${myArray[0]}_simulation_$p.fastq
  done
done < "$tinput"
echo "....Completed mixtures."
#note: mixtures are named by the subject


# --- MetaSPAdes & bowtie2-build per donor / pre file (8.75M depth) ---
echo "--- Running metaSPAdes and bowtie2-build for databases...."
module load bowtie2
# Loop over subjects in triads table
tinput="./triads_table.txt"
i=1
while ICF=$'\t' read -r -a myArray
do
  # skip header for now
  test $i -eq 1 && ((i=i+1)) && continue
  # metaSPAdes per donor and pre kneaded FASTQs (full depth)
  echo "Running metaSPAdes for ${myArray[0]} : donor_${myArray[3]} pre_${myArray[2]}"
  ## donor file
  if [ ! -e "metaspades_out/donor_${myArray[3]}" ]; then
    /home/knightsd/shoops/lib/SPAdes-3.15.5-Linux/bin/metaspades.py -t 128 -1 ./PRJEB23524/kneaded/${myArray[9]}_kneaddata_paired_1.fastq -2 ./PRJEB23524/kneaded/${myArray[9]}_kneaddata_paired_2.fastq -o ./metaspades_out/donor_${myArray[3]}
  fi
  ## pre file
  /home/knightsd/shoops/lib/SPAdes-3.15.5-Linux/bin/metaspades.py -t 128 -1 ./PRJEB23524/kneaded/${myArray[8]}_kneaddata_paired_1.fastq -2 ./PRJEB23524/kneaded/${myArray[8]}_kneaddata_paired_2.fastq -o ./metaspades_out/pre_${myArray[2]}
  # bowtie2-build from output
  echo "Building databases for ${myArray[0]} : donor_${myArray[3]} pre_${myArray[2]}"
  ## donor MAGs
  bowtie2-build --seed 25 --threads 16 metaspades_out/donor_${myArray[3]}/scaffolds.fasta bowtie_db/donor_${myArray[3]}
  ## pre MAGs
  bowtie2-build --seed 25 --threads 16 metaspades_out/pre_${myArray[2]}/scaffolds.fasta bowtie_db/pre_${myArray[2]}
done < "$tinput"
echo "....Completed databases."


# --- Align Mixtures with databases ---
echo "--- Aligning of mixtures and databases...."
# Loop over each donor:pre pair
tinput="./triads_table.txt"
i=1
while ICF=$'\t' read -r -a myArray
do
  # skip header for now
  test $i -eq 1 && ((i=i+1)) && continue
  # Loop over simulation amounts
  for p in 0 10 20 30 40 50 60 70 80 90 100
  do
    # 1.00 identity alignments
    /usr/bin/time bowtie2-align-s -x ./bowtie_db/donor_${myArray[3]} -S ./bowtie_out/${myArray[0]}_sim_$p_donor_alignment_1.00.txt --very-sensitive --np 1 --mp "1,1" --rdg "0,1" --rfg "0,1" --score-min "C,0" -k 1  -p 16 --no-hd --no-unal -q ./mixtures/${myArray[0]}_simulation_$p.fastq
    /usr/bin/time bowtie2-align-s -x ./bowtie_db/pre_${myArray[2]} -S ./bowtie_out/${myArray[0]}_sim_$p_pre_alignment_1.00.txt --very-sensitive --np 1 --mp "1,1" --rdg "0,1" --rfg "0,1" --score-min "C,0" -k 1  -p 16 --no-hd --no-unal -q ./mixtures/${myArray[0]}_simulation_$p.fastq
    # 0.99 identity alignments
    /usr/bin/time bowtie2-align-s -x ./bowtie_db/donor_${myArray[3]} -S ./bowtie_out/${myArray[0]}_sim_$p_donor_alignment_0.99.txt --very-sensitive --np 1 --mp "1,1" --rdg "0,1" --rfg "0,1" --score-min "L,0,-.01" -k 1  -p 16 --no-hd --no-unal -q ./mixtures/${myArray[0]}_simulation_$p.fastq
    /usr/bin/time bowtie2-align-s -x ./bowtie_db/pre_${myArray[2]} -S ./bowtie_out/${myArray[0]}_sim_$p_pre_alignment_0.99.txt --very-sensitive --np 1 --mp "1,1" --rdg "0,1" --rfg "0,1" --score-min "L,0,-.01" -k 1  -p 16 --no-hd --no-unal -q ./mixtures/${myArray[0]}_simulation_$p.fastq
    # 0.98 identity alignments
    /usr/bin/time bowtie2-align-s -x ./bowtie_db/donor_${myArray[3]} -S ./bowtie_out/${myArray[0]}_sim_$p_donor_alignment_0.98.txt --very-sensitive --np 1 --mp "1,1" --rdg "0,1" --rfg "0,1" --score-min "L,0,-.02" -k 1  -p 16 --no-hd --no-unal -q ./mixtures/${myArray[0]}_simulation_$p.fastq
    /usr/bin/time bowtie2-align-s -x ./bowtie_db/pre_${myArray[2]} -S ./bowtie_out/${myArray[0]}_sim_$p_pre_alignment_0.98.txt --very-sensitive --np 1 --mp "1,1" --rdg "0,1" --rfg "0,1" --score-min "L,0,-.02" -k 1  -p 16 --no-hd --no-unal -q ./mixtures/${myArray[0]}_simulation_$p.fastq
  done
done < "$tinput"
echo "....Completed alignment."


# --- Mixtures & alignments for 20% contamination ---
echo "--- Contaminate mixtures & alignments...."
# Loop over each donor:pre pair
tinput="./triads_table.txt"
i=1
while ICF=$'\t' read -r -a myArray
do
  # skip header for now
  test $i -eq 1 && ((i=i+1)) && continue
  # Loop over mixtures
  for p in 0 10 20 30 40 50 60 70 80 90 100
  do
    # Create Mixture
    ## initialize with 20% contaminate
    head -n 1750000 PRJEB23524/kneaded/${myArray[13]}_kneaddata_paired_1.fastq > ./contaminated/${myArray[0]}_simulation_$p_extra20.fastq
    ## add in remaining donor
    donor_share=`expr 70000 \\* $p`
    head -n $donor_share PRJEB23524/kneaded/${myArray[9]}_kneaddata_paired_1.fastq >> ./contaminated/${myArray[0]}_simulation_$p_extra20.fastq
    ## add in remaining pre
    invp=`expr 100 - $p`
    pre_share=`expr 70000 \\* $invp`
    head -n $pre_share PRJEB23524/kneaded/${myArray[8]}_kneaddata_paired_1.fastq >> ./contaminated/${myArray[0]}_simulation_$p_extra20.fastq
    # Align mixture with respective databases
    ## 1.00 identity alignments
    /usr/bin/time bowtie2-align-s -x ./bowtie_db/donor_${myArray[3]} -S ./contaminated/bowtie_out/${myArray[0]}_sim_$p_extra20_donor_alignment_1.00.txt --very-sensitive --np 1 --mp "1,1" --rdg "0,1" --rfg "0,1" --score-min "C,0" -k 1  -p 16 --no-hd --no-unal -q ./contaminated/${myArray[0]}_simulation_$p_extra20.fastq
    /usr/bin/time bowtie2-align-s -x ./bowtie_db/pre_${myArray[2]} -S ./contaminated/bowtie_out/${myArray[0]}_sim_$p_extra20_pre_alignment_1.00.txt --very-sensitive --np 1 --mp "1,1" --rdg "0,1" --rfg "0,1" --score-min "C,0" -k 1  -p 16 --no-hd --no-unal -q ./contaminated/${myArray[0]}_simulation_$p_extra20.fastq
    ## 0.99 identity alignments
    /usr/bin/time bowtie2-align-s -x ./bowtie_db/donor_${myArray[3]} -S ./contaminated/bowtie_out/${myArray[0]}_sim_$p_extra20_donor_alignment_0.99.txt --very-sensitive --np 1 --mp "1,1" --rdg "0,1" --rfg "0,1" --score-min "L,0,-.01" -k 1  -p 16 --no-hd --no-unal -q ./contaminated/${myArray[0]}_simulation_$p_extra20.fastq
    /usr/bin/time bowtie2-align-s -x ./bowtie_db/pre_${myArray[2]} -S ./contaminated/bowtie_out/${myArray[0]}_sim_$p_extra20_pre_alignment_0.99.txt --very-sensitive --np 1 --mp "1,1" --rdg "0,1" --rfg "0,1" --score-min "L,0,-.01" -k 1  -p 16 --no-hd --no-unal -q ./contaminated/${myArray[0]}_simulation_$p_extra20.fastq
    ## 0.98 identity alignments
    /usr/bin/time bowtie2-align-s -x ./bowtie_db/donor_${myArray[3]} -S ./contaminated/bowtie_out/${myArray[0]}_sim_$p_extra20_donor_alignment_0.98.txt --very-sensitive --np 1 --mp "1,1" --rdg "0,1" --rfg "0,1" --score-min "L,0,-.02" -k 1  -p 16 --no-hd --no-unal -q ./contaminated/${myArray[0]}_simulation_$p_extra20.fastq
    /usr/bin/time bowtie2-align-s -x ./bowtie_db/pre_${myArray[2]} -S ./contaminated/bowtie_out/${myArray[0]}_sim_$p_extra20_pre_alignment_0.98.txt --very-sensitive --np 1 --mp "1,1" --rdg "0,1" --rfg "0,1" --score-min "L,0,-.02" -k 1  -p 16 --no-hd --no-unal -q ./contaminated/${myArray[0]}_simulation_$p_extra20.fastq
  done
done < "$tinput"
echo "....Completed extra20 mixtures and alignments."


# --- Parse output per pair ---
echo "--- Parsing output files...."
tinput="./triads_table.txt"
i=1
while ICF=$'\t' read -r -a myArray
do
  # skip header for now
  test $i -eq 1 && ((i=i+1)) && continue
  # per mixture
  for p in 0 10 20 30 40 50 60 70 80 90 100
  do
    # Regular Mixtures
    ## 1.00 identity counts
    python3 parse_sim_output.py ${myArray[0]} $p 1.00 counts_simulations_1.00.txt
    ## 0.99 identity counts
    python3 parse_sim_output.py ${myArray[0]} $p 0.99 counts_simulations_0.99.txt
    ## 0.98 identity counts
    python3 parse_sim_output.py ${myArray[0]} $p 0.98 counts_simulations_0.98.txt
    # Contaminated mixtures
    ## 1.00 identity counts
    python3 parse_contam_output.py ${myArray[0]} $p 1.00 counts_contam_sim_1.00.txt
    ## 0.99 identity counts
    python3 parse_contam_output.py ${myArray[0]} $p 0.99 counts_contam_sim_0.99.txt
    ## 0.98 identity counts
    python3 parse_contam_output.py ${myArray[0]} $p 0.98 counts_contam_sim_0.98.txt
  done
done < "$tinput"
echo "....Completed parsing."

# --- Conclusion ---
echo "All done."
