#!/bin/bash
#SBATCH -J AMR++ -o AMR++_log.out -t 24:00:00 -p interactive --mem=200G --nodes=1 --ntasks=4 --cpus-per-task=8

# Remember to change the "maxForks" parameter in the config file that you are using. This corresponds with "--ntasks"
# to control how many jobs are submitted at once. The "--threads" argument should also match the --cpus-per-task to
# fully utlize the available computing resources.

# This script works on TAMU's Grace HPRC, but you need to follow the instructions on the Github to get the right conda
# environment installed on your computing environment
source ~/.bashrc
conda activate AMR++_env  # Explore the installation instructions on github to see how to install this environment

# Because we are using the local profile, software from the conda environment will be used. The flag "maxForks" is
# in the local.config file and will spawn 4 processes by default, this corresponds with "--ntasks=4" in the sbatch script.
#nextflow run main_AMR++.nf -profile local --threads 8 # This will use 8 threads, which corresponds with "--cpus-per-task=8".

# SUZIE ADDED THESE LINES TO FOLLOW
cd ~/downloads/AMRplusplus/
# Evaluation QC pipeline
#nextflow run main_AMR++.nf -profile local --threads 8 \
#  --pipeline eval_qc \
#  --reads "/scratch.global/shoops/magic_gbs_068/*R{1,2}_001.fastq.gz" \
#  --output "/scratch.global/shoops/amr_out/"
#echo "Completed QC evaluation."

# Trimming QC pipeline
#nextflow run main_AMR++.nf -profile local --threads 8 \
#  --pipeline trim_qc \
#  --reads "/scratch.global/shoops/magic_gbs_068/*R{1,2}_001.fastq.gz" \
#  --output "/scratch.global/shoops/amr_out/"
#echo "Completed trimming."

# Host removal QC pipeline
nextflow run main_AMR++.nf -profile local --threads 8 \
  --pipeline rm_host \
  --reads "/scratch.global/shoops/amr_out/QC_trimming/Paired/*{1,2}P.fastq.gz" \
  --output "/scratch.global/shoops/amr_out/"
echo "Completed host removal."

# Resistome pipeline
nextflow run main_AMR++.nf -profile local --threads 8 \
  --pipeline resistome \
  --reads "/scratch.global/shoops/amr_out/HostRemoval/NonHostFastq/*R{1,2}.fastq.gz" \
  --amr "~/downloads/AMRplusplus/data/amr/megares_database_v3.00.fasta" \
  --annotation "~/downloads/AMRplusplus/data/amr/megares_annotations_v3.00.csv" \
  --output "/scratch.global/shoops/amr_out/"
echo "Completed resistome alignment."

echo "All done."