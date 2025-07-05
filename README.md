# MAGEnTa

## Installation instructions
Install necessary software and python packages
### Recommended : New conda environment
If you use conda, we suggest a new virtual environment using python version 3.7
```
conda create --name magenta python=3.7
```

### Python Packages
Install the following python packages. We recommend using pip to install the packages within your virtual environment.

Required packages include: os (included in python), sys (included in python), numpy
```
pip install numpy
```

### Software Installation
metaSPAdes : Follow metaspades instructions outlined here: https://ablab.github.io/spades/installation.html#downloading-spades-linux-binaries
The SPAdes library should be automatically added to your PATH after you compile.
Bowtie2 : install bowtie2 using instruction here (Conda method suggested): https://www.metagenomics.wiki/tools/bowtie2/install
add here: check if these softwares are available/installed


## Usage and Tutorial
### Set-Up
See the above installation instructions for setting up your virtual environment and installing all the necessary dependencies. The instructions in this tutorial assume you have cloned this repo and followed installation instructions.
```
# Go the the cloned repo path
cd /path/to/magenta/repo
# Activate your virtual environment (optional)
conda activate magenta
```

### Prepare Sample Table
The MAGEnTa tool needs a **table containing the sample identifiers in a pre-defined format**. For this tutorial, we have a provided an example table (**sim_magenta_table.txt**), which users may also use a template for preparing a table for their own dataset.

The tutorial sample table is shown below:

** | Post_ID | Post_Filename | Donor_ID | Donor_Filename_fwd | Donor_Filename_rev | Pre_ID | Pre_Filename_fwd | Pre_Filename_rev | **
| ------- | ------------- | -------- | ------------------ | ------------------ | ------ | ---------------- | ---------------- |
| Sim0Pre | simulation_0.0_pre.fastq | Donor1 | donor_simulation_r1.fastq | donor_simulation_r2.fastq | Pre1 | pre_simulation_r1.fastq | pre_simulation_r2.fastq |
| Sim10Pre | simulation_0.1_pre.fastq | Donor1 | donor_simulation_r1.fastq | donor_simulation_r2.fastq | Pre1 | pre_simulation_r1.fastq | pre_simulation_r2.fastq |
| Sim20Pre | simulation_0.2_pre.fastq | Donor1 | donor_simulation_r1.fastq | donor_simulation_r2.fastq | Pre1 | pre_simulation_r1.fastq | pre_simulation_r2.fastq |
| Sim30Pre | simulation_0.3_pre.fastq | Donor1 | donor_simulation_r1.fastq | donor_simulation_r2.fastq | Pre1 | pre_simulation_r1.fastq | pre_simulation_r2.fastq |
| Sim40Pre | simulation_0.4_pre.fastq | Donor1 | donor_simulation_r1.fastq | donor_simulation_r2.fastq | Pre1 | pre_simulation_r1.fastq | pre_simulation_r2.fastq |
| Sim50Pre | simulation_0.5_pre.fastq | Donor1 | donor_simulation_r1.fastq | donor_simulation_r2.fastq | Pre1 | pre_simulation_r1.fastq | pre_simulation_r2.fastq |
| Sim60Pre | simulation_0.6_pre.fastq | Donor1 | donor_simulation_r1.fastq | donor_simulation_r2.fastq | Pre1 | pre_simulation_r1.fastq | pre_simulation_r2.fastq |
| Sim70Pre | simulation_0.7_pre.fastq | Donor1 | donor_simulation_r1.fastq | donor_simulation_r2.fastq | Pre1 | pre_simulation_r1.fastq | pre_simulation_r2.fastq |
| Sim80Pre | simulation_0.8_pre.fastq | Donor1 | donor_simulation_r1.fastq | donor_simulation_r2.fastq | Pre1 | pre_simulation_r1.fastq | pre_simulation_r2.fastq |
| Sim90Pre | simulation_0.9_pre.fastq | Donor1 | donor_simulation_r1.fastq | donor_simulation_r2.fastq | Pre1 | pre_simulation_r1.fastq | pre_simulation_r2.fastq |
| Sim100Pre | simulation_1.0_pre.fastq | Donor1 | donor_simulation_r1.fastq | donor_simulation_r2.fastq | Pre1 | pre_simulation_r1.fastq | pre_simulation_r2.fastq |


### MAGs Database Build
Produces MAGs databases per donor and pre-treatment sample file.
```
# Produce MAGs databases (most time-consuming step, may need to run on multiple threads and use a cluster manager)
./src/build_mag_db.sh -i INPUT_DIR -m MAPPING_FILE -o OUTPUT_DIR -t THREADS
```

Instructions for the tutorial:
```
./src/build_mag_db.sh -i data/sim_samples -m data/sim_magenta_table.txt -o data/sim_results -t 256
```


### Alignmeant with MAGs Databases
Produces counts table and engrafted mags lists.
```
# Align the samples with
./src/align_source.sh -i INPUT_DIR -m MAPPING_FILE -o OUTPUT_DIR -t THREADS -s ALIGN_SCORE
```
note that "ALIGN_SCORE" can only be one of 3 choices for alignment score in the current implementation: 1.00, 0.99, 0.98

Instructions for the tutorial:
```
./src/align_source.sh -i sim_results -m data/sim_magenta_table.txt -o data/sim_results -t 16 -s 1.00
```

### Tabulate the alignment output to get Engraftment Estimates
The python script **parse_bt2_output.py** evaluates the alignment files to provide engraftement estimates per sample set.
```
# Get engraftment tables output
python3 ./src/parse_bt2_output.py -i [alignment_output_path] -m [mapping_file] -o [output_table]
```

To follow the tutorial:
```
parse_bt2_output.py -i data/sim_results/align -m data/sim_magenta_table.txt -o sim_results/sim_out_table.txt
```


