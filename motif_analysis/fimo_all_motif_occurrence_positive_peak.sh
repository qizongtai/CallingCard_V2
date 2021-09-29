#!/bin/bash
#
#SBATCH --job-name=get_motif_per_peak
#SBATCH --output=logs/get_motif_per_peak_%A_%a.out
#SBATCH --error=logs/get_motif_per_peak_%A_%a.err
#SBATCH --array=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=16000

###---get_motif_occurrence_per_peak---###
tf[1]="SNAI2"
tf[2]="FOSL1"
tf[3]="FOS_JUN"

cell[1]="SCC9"
cell[2]="JHU6"
cell[3]="SCC47"
cell[4]="SCC4"

out_dir="../3_output/fimo_all_motif_occurance_peak_all_a"
[ ! -d $out_dir ] && mkdir $out_dir
out_file="$out_dir/summary_peak_with_motif_occurrence.txt"
[ -e $out_file ] && rm $out_file

for i in {1..3}; do
    out_dir_1="../3_output/fimo_${tf[$i]}_SNAI2"
    [ ! -d $out_dir_1 ] && mkdir $out_dir_1
    out_file_1="$out_dir_1/summary_peak_with_motif_occurrence.txt"
    [ -e $out_file_1 ] && rm $out_file_1
    out_file_2="$out_dir_1/summary_peak_with_motif_occurrence_extracted.txt"
    [ -e $out_file_2 ] && rm $out_file_2
    out_file_3="$out_dir_1/summary_peak_with_motif_occurrence_rowname.txt"
    [ -e $out_file_3 ] && rm $out_file_3
    out_file_4="$out_dir_1/summary_peak_with_motif_occurrence_extracted_rowname.txt"
    [ -e $out_file_4 ] && rm $out_file_4
    for k in {1..4}; do
        #create occurrence_per_peak.txt file for each motif, each scale and each cell
        out_file_a1="../3_output/fimo_${tf[$i]}_SNAI2_${cell[$k]}_macs_win500_step250_p0001/fimo_${tf[$i]}_SNAI2_${cell[$k]}_occurrence_per_peak.txt"
        [ -e $out_file_a1 ] && rm $out_file_a1
        in_file_a1="../3_output/fimo_${tf[$i]}_SNAI2_${cell[$k]}_macs_win500_step250_p0001/fimo.tsv"
        cut -f 3 $in_file_a1 | sort | uniq -c | sed -r 's/\s+/\t/g' | sed -r 's/^\t//' | tail -n +2 | sort -nr -k 1 | egrep -v "(#|sequence_name)" > $out_file_a1
        #extract and add to the summary file
        wc -l $out_file_a1 >> $out_file_1
    done
    #creat rowname
    echo "${tf[$i]}" >> $out_file_3

    #extracted and formatted
    cat $out_file_1 | cut -f 1 -d ' ' | xargs -n 4 | sed 's/ /\t/g' >  $out_file_2
    #merge rowname with occurrence_per_peak_extracted.txt
    paste $out_file_3 $out_file_2 > $out_file_4
    #cat all tf data to one file by >>
    paste $out_file_3 $out_file_2 >> $out_file
done


