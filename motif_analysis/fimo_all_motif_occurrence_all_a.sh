#!/bin/bash
#
#SBATCH --job-name=bash
#SBATCH --output=logs/bash_%A_%a.out
#SBATCH --error=logs/bash_%A_%a.err
#SBATCH --array=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=16000

###---data extraction---###
#:<<'END'
tf[1]="SNAI2"
tf[2]="SNAI1"
tf[3]="ZEB1"
tf[4]="TWIST1"
tf[5]="TWIST2"
tf[6]="FOSL1"
tf[7]="BATF"
tf[8]="BATF3"
tf[9]="FOS_JUN"
tf[10]="BATF_JUN"
tf[11]="CTCF"
tf[12]="OLIG2"
tf[13]="SP1"
tf[14]="SOX10"
tf[15]="HINFP"

scale[1]="05"
scale[2]="01"
scale[3]="005"
scale[4]="001"
scale[5]="0005"
scale[6]="0001"
scale[7]="00005"
scale[8]="00001"
scale[9]="000005"
scale[10]="000001"
scale[11]="0000005"
scale[12]="0000001"

out_dir="../3_output/fimo_all_motif_occurrence_all_a"
out_file_summary="$out_dir/summary_extracted_rowname.txt"
[ ! -d $out_dir ] && mkdir $out_dir
[ -e $out_file_summary ] && rm $out_file_summary

for i in {1..15}; do
     out_file_rowname="../3_output/fimo_${tf[$i]}_SNAI2/rowname.txt"
     [ -e $out_file_rowname ] && rm $out_file_rowname
     for j in {1..12}; do
          #create rowname file at each a and tf
          echo "${tf[$i]}_${scale[$j]}" >> $out_file_rowname
     done
     #locate the data directory and add the row names
     #paste $out_file_rowname "../3_output/fimo_${tf[$i]}_SNAI2/summary_extracted.txt" > ../3_output/fimo_${tf[$i]}_SNAI2/summary_extracted_rowname.txt
     paste $out_file_rowname "../3_output/fimo_${tf[$i]}_SNAI2/summary_extracted.txt" >> $out_file_summary
done

