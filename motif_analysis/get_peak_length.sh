#script extracted from
#/scratch/splab/zqi/callingcards_bulk/7_BulkCC_hyPBase_SNAI2_SCC47_MGI/2_peak_caller/2_code/1_analysis_202007.sh

###--------------extracted----------------###
#!/bin/bash
#
#SBATCH --job-name=analysis
#SBATCH --output=logs/analysis_%A_%a.out
#SBATCH --error=logs/analysis_%A_%a.err
#SBATCH --array=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=1000

###---calculate entire length---###
#:<<'END'
in_dir="/scratch/splab/zqi/callingcard_bulk/10_BulkCC_SNAI2_different_peak_calling/3_output/macs_win500_step250/macs"

in_file[1]="${in_dir}/SNAI2_SCC9_clean_peak_macs_win500_step250_p0001.peaks.bed"
in_file[2]="${in_dir}/SNAI2_JHU6_clean_peak_macs_win500_step250_p0001.peaks.bed"
in_file[3]="${in_dir}/SNAI2_SCC47_clean_peak_macs_win500_step250_p0001.peaks.bed"
in_file[4]="${in_dir}/SNAI2_SCC4_clean_peak_macs_win500_step250_p0001.peaks.bed"


out_dir="../3_output/"
out_file="SNAI2_peak_macs_win500_step250_p0001_length_summary.txt"
[ -e ${out_dir}/${out_file} ] && rm ${out_dir}/${out_file}
for i in {1..4}; do
  awk 'BEGIN{OFS=FS="\t"} {sum += ($3-$2)} END{print sum}' ${in_file[$i]} >> ${out_dir}/${out_file}
  done

:<<'END'
#make a 4-column format
cat ${out_dir}/${out_file} | xargs -n 4 | sed 's/ /\t/g' >  $out_dir/SNAI2_peak_length_summary_formatted.txt
#make 15 copies vertically in one file
copy=$out_dir/SNAI2_peak_length_summary_formatted.txt
cat $copy $copy $copy $copy $copy $copy $copy $copy $copy $copy $copy $copy $copy $copy $copy > $out_dir/SNAI2_peak_length_summary_formatted_15copies.txt
END
