#!/bin/bash
#
#SBATCH --job-name=get_peak_number
#SBATCH --output=logs/get_peak_number_%A_%a.out
#SBATCH --error=logs/get_peak_number_%A_%a.err
#SBATCH --array=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=16000

###---get peak number---###
in_file[1]="../3_output/SNAI2_SCC9_a05_bonf_clean.fa"
in_file[2]="../3_output/SNAI2_JHU6_a05_bonf_clean.fa"
in_file[3]="../3_output/SNAI2_SCC47_a05_bonf_clean.fa"
in_file[4]="../3_output/SNAI2_SCC4_a05_bonf_clean.fa"
in_file[5]="../3_output/SNAI2_SCC9_a01_bonf_clean.fa"
in_file[6]="../3_output/SNAI2_JHU6_a01_bonf_clean.fa"
in_file[7]="../3_output/SNAI2_SCC47_a01_bonf_clean.fa"
in_file[8]="../3_output/SNAI2_SCC4_a01_bonf_clean.fa"
in_file[9]="../3_output/SNAI2_SCC9_a005_bonf_clean.fa"
in_file[10]="../3_output/SNAI2_JHU6_a005_bonf_clean.fa"
in_file[11]="../3_output/SNAI2_SCC47_a005_bonf_clean.fa"
in_file[12]="../3_output/SNAI2_SCC4_a005_bonf_clean.fa"
in_file[13]="../3_output/SNAI2_SCC9_a001_bonf_clean.fa"
in_file[14]="../3_output/SNAI2_JHU6_a001_bonf_clean.fa"
in_file[15]="../3_output/SNAI2_SCC47_a001_bonf_clean.fa"
in_file[16]="../3_output/SNAI2_SCC4_a001_bonf_clean.fa"
in_file[17]="../3_output/SNAI2_SCC9_a0005_bonf_clean.fa"
in_file[18]="../3_output/SNAI2_JHU6_a0005_bonf_clean.fa"
in_file[19]="../3_output/SNAI2_SCC47_a0005_bonf_clean.fa"
in_file[20]="../3_output/SNAI2_SCC4_a0005_bonf_clean.fa"
in_file[21]="../3_output/SNAI2_SCC9_a0001_bonf_clean.fa"
in_file[22]="../3_output/SNAI2_JHU6_a0001_bonf_clean.fa"
in_file[23]="../3_output/SNAI2_SCC47_a0001_bonf_clean.fa"
in_file[24]="../3_output/SNAI2_SCC4_a0001_bonf_clean.fa"
in_file[25]="../3_output/SNAI2_SCC9_a00005_bonf_clean.fa"
in_file[26]="../3_output/SNAI2_JHU6_a00005_bonf_clean.fa"
in_file[27]="../3_output/SNAI2_SCC47_a00005_bonf_clean.fa"
in_file[28]="../3_output/SNAI2_SCC4_a00005_bonf_clean.fa"
in_file[29]="../3_output/SNAI2_SCC9_a00001_bonf_clean.fa"
in_file[30]="../3_output/SNAI2_JHU6_a00001_bonf_clean.fa"
in_file[31]="../3_output/SNAI2_SCC47_a00001_bonf_clean.fa"
in_file[32]="../3_output/SNAI2_SCC4_a00001_bonf_clean.fa"
in_file[33]="../3_output/SNAI2_SCC9_a000005_bonf_clean.fa"
in_file[34]="../3_output/SNAI2_JHU6_a000005_bonf_clean.fa"
in_file[35]="../3_output/SNAI2_SCC47_a000005_bonf_clean.fa"
in_file[36]="../3_output/SNAI2_SCC4_a000005_bonf_clean.fa"
in_file[37]="../3_output/SNAI2_SCC9_a000001_bonf_clean.fa"
in_file[38]="../3_output/SNAI2_JHU6_a000001_bonf_clean.fa"
in_file[39]="../3_output/SNAI2_SCC47_a000001_bonf_clean.fa"
in_file[40]="../3_output/SNAI2_SCC4_a000001_bonf_clean.fa"
in_file[41]="../3_output/SNAI2_SCC9_a0000005_bonf_clean.fa"
in_file[42]="../3_output/SNAI2_JHU6_a0000005_bonf_clean.fa"
in_file[43]="../3_output/SNAI2_SCC47_a0000005_bonf_clean.fa"
in_file[44]="../3_output/SNAI2_SCC4_a0000005_bonf_clean.fa"
in_file[45]="../3_output/SNAI2_SCC9_a0000001_bonf_clean.fa"
in_file[46]="../3_output/SNAI2_JHU6_a0000001_bonf_clean.fa"
in_file[47]="../3_output/SNAI2_SCC47_a0000001_bonf_clean.fa"
in_file[48]="../3_output/SNAI2_SCC4_a0000001_bonf_clean.fa"

#option1: output is NOT ordered as they run in parallel and may ends at different times
#id=${SLURM_ARRAY_TASK_ID}

#option2: output is ordered
out_sum_dir="../3_output/peak_SNAI2_all_a"
[ ! -d $out_sum_dir ] && mkdir $out_sum_dir
for id in {1..48}; do
     wc -l ${in_file[$id]} >> $out_sum_dir/summary.txt
     done
#extract numbers and make 4 column table
#cat $out_sum_dir/summary.txt | cut -f 1 -d ' ' | xargs -n 4 | sed 's/ /\t/g' >  $out_sum_dir/summary_extracted_2.txt
cat $out_sum_dir/summary.txt | cut -f 1 -d ' ' | awk '{print $1/2}' | xargs -n 4 | sed 's/ /\t/g' >  $out_sum_dir/summary_extracted.txt
copy="$out_sum_dir/summary_extracted.txt"
cat $copy $copy $copy $copy $copy $copy $copy $copy $copy $copy $copy $copy $copy $copy $copy > $out_sum_dir/summary_extracted_15copies.txt
#END