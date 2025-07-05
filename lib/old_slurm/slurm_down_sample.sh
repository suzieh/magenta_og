##### Down Sampling #####
# Using seqtk to down sample fastqs (be sure to set seed for same reads in each direction):
./seqtk/seqtk sample -s125 input_R1.fastq 500000 > down_500k_R2.fastq

# Have a copy of kneaddata output in ucfmt_kneaddata
cp kneaddata_out/*paired_[1-2].fastq ucfmt_kneaddata

# Have a temporary holder for all kneaddata output
cp ucfmt_kneaddata/prep_samples/*paired_[1-2].fastq kneaddata_out

# Loop(s) to downsample to each level
#     note: kneaddata output is not in the same order, so we have to use subseq in order to match ids
# 500k attempt
cd kneaddata_out/
for fn in *_paired_1.fastq
do
  ../seqtk/seqtk sample -s125 $fn 500000 > ../ucfmt_500k/$fn
  grep @A00223 ../ucfmt_500k/$fn > tmp_ids.txt
  sed -i 's/\/1/\/2/g' tmp_ids.txt
  sed -i 's/^@//g' tmp_ids.txt
  ../seqtk/seqtk subseq ${fn/paired_1/paired_2} tmp_ids.txt > ../ucfmt_500k/${fn/paired_1/paired_2}
done

# 1M attempt
cd kneaddata_out/
for fn in *_paired_1.fastq
do
  ../seqtk/seqtk sample -s125 $fn 1000000 > ../ucfmt_1M/$fn
  grep @A00223 ../ucfmt_1M/$fn > tmp_ids.txt
  sed -i 's/\/1/\/2/g' tmp_ids.txt
  sed -i 's/^@//g' tmp_ids.txt
  ../seqtk/seqtk subseq ${fn/paired_1/paired_2} tmp_ids.txt > ../ucfmt_1M/${fn/paired_1/paired_2}
done

# 2M attempt
cd kneaddata_out/
for fn in *_paired_1.fastq
do
  ../seqtk/seqtk sample -s125 $fn 2000000 > ../ucfmt_2M/$fn
  grep @A00223 ../ucfmt_2M/$fn > tmp_ids_2M.txt
  sed -i 's/\/1/\/2/g' tmp_ids_2M.txt
  sed -i 's/^@//g' tmp_ids_2M.txt
  ../seqtk/seqtk subseq ${fn/paired_1/paired_2} tmp_ids_2M.txt > ../ucfmt_2M/${fn/paired_1/paired_2}
done

# 5M attempt
cd kneaddata_out/
for fn in *_paired_1.fastq
do
  ../seqtk/seqtk sample -s125 $fn 5000000 > ../ucfmt_5M/$fn
  grep @A00223 ../ucfmt_5M/$fn > tmp_ids_5M.txt
  sed -i 's/\/1/\/2/g' tmp_ids_5M.txt
  sed -i 's/^@//g' tmp_ids_5M.txt
  ../seqtk/seqtk subseq ${fn/paired_1/paired_2} tmp_ids_5M.txt > ../ucfmt_5M/${fn/paired_1/paired_2}
done

# 10M attempt
cd kneaddata_out/
for fn in *_paired_1.fastq
do
  ../seqtk/seqtk sample -s125 $fn 10000000 > ../ucfmt_10M/$fn
  grep @A00223 ../ucfmt_10M/$fn > tmp_ids_10M.txt
  sed -i 's/\/1/\/2/g' tmp_ids_10M.txt
  sed -i 's/^@//g' tmp_ids_10M.txt
  ../seqtk/seqtk subseq ${fn/paired_1/paired_2} tmp_ids_10M.txt > ../ucfmt_10M/${fn/paired_1/paired_2}
done


# TEMP: testing the approach to get the same subset of samples
# cd test_down
# ../seqtk/seqtk sample -s125 UCFMT011BL_S66_R1_001_kneaddata_paired_1.fastq 500000 > down_500k_UCFMT011BL_S66_R1_001_kneaddata_paired_1.fastq
# grep @A00223 down_500k_UCFMT011BL_S66_R1_001_kneaddata_paired_1.fastq > tmp_ids.txt
# sed -i 's/\/1/\/2/g' tmp_ids.txt
# sed -i 's/^@//g' tmp_ids.txt
# ../seqtk/seqtk subseq UCFMT011BL_S66_R1_001_kneaddata_paired_2.fastq tmp_ids.txt > down_500k_UCFMT011BL_S66_R1_001_kneaddata_paired_2.fastq




