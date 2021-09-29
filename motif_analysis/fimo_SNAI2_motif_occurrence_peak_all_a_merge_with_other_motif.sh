#!/bin/bash
#
#SBATCH --job-name=motif_merge
#SBATCH --output=logs/motif_merge_%A_%a.out
#SBATCH --error=logs/motif_merge_%A_%a.err
#SBATCH --array=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=10G

###---get_motif_occurrence_per_peak---###
tf[1]="SNAI1"
tf[2]="ZEB1"
tf[3]="TWIST1"
tf[4]="TWIST2"
tf[5]="FOSL1"
tf[6]="BATF"
tf[7]="BATF3"
tf[8]="FOS_JUN"
tf[9]="BATF_JUN"
tf[10]="CTCF"
tf[11]="OLIG2"
tf[12]="SP1"
tf[13]="SOX10"
tf[14]="HINFP"

scale[1]="a05"
scale[2]="a01"
scale[3]="a005"
scale[4]="a001"
scale[5]="a0005"
scale[6]="a0001"
scale[7]="a00005"
scale[8]="a00001"
scale[9]="a000005"
scale[10]="a000001"
scale[11]="a0000005"
scale[12]="a0000001"

cell[1]="SCC9"
cell[2]="JHU6"
cell[3]="SCC47"
cell[4]="SCC4"

out_dir="../3_output/fimo_SNAI2_motif_occurrence_peak_all_a_merge_with_other_motif"
[ ! -d $out_dir ] && mkdir $out_dir
out_file1="$out_dir/fimo_SNAI2_SNAI2_vs_other_motif_summary.txt"
[ -e $out_file1 ] && rm $out_file1
out_file2="$out_dir/fimo_SNAI2_SNAI2_vs_other_motif_rowname_all_level.txt"
[ -e $out_file2 ] && rm $out_file2
out_file3="$out_dir/fimo_SNAI2_SNAI2_vs_other_motif_rowname_collapse_cell.txt"
[ -e $out_file3 ] && rm $out_file3

for i in {1..14}; do
    for j in {1..12}; do
        for k in {1..4}; do
            #create occurrence_per_peak.txt file for each motif, each scale and each cell
            out_file6="$out_dir/fimo_SNAI2_SNAI2_vs_${tf[$i]}_${cell[$k]}_${scale[$j]}.txt"
            out_file7="$out_dir/fimo_SNAI2_SNAI2_vs_${tf[$i]}_${cell[$k]}_${scale[$j]}.bed"
            in_file_1="../3_output/fimo_SNAI2_SNAI2_${cell[$k]}_${scale[$j]}/fimo_SNAI2_SNAI2_${cell[$k]}_${scale[$j]}_postive_peak_and_occurrence_per_peak.txt"
            in_file_2="../3_output/fimo_${tf[$i]}_SNAI2_${cell[$k]}_${scale[$j]}/fimo_${tf[$i]}_SNAI2_${cell[$k]}_${scale[$j]}_occurrence_per_peak.txt"
            join -1 1 -2 2 -t $'\t' <(sort -k 1 $in_file_1) <(sort -k 2 $in_file_2) > $out_file6
            awk 'BEGIN{OFS=FS="\t"} {print $3,$4,$5,$6,$7,$8,$1,$2,$9}' $out_file6 | sort -n -r -k 8 > $out_file7
            #extract and add to the summary file
            wc -l $in_file_1 >> $out_file1
            wc -l $in_file_2 >> $out_file1
            wc -l $out_file6 >> $out_file1
            #creat rowname at all levels
            echo "SANI2_vs_${tf[$i]}_${cell[$k]}_${scale[$j]}" >> $out_file2
        done
        #creat rowname collapse on cells
        echo "SANI2_vs_${tf[$i]}_${scale[$j]}" >> $out_file3
    done
done

out_file8="$out_dir/fimo_SNAI2_SNAI2_vs_other_motif_summary_all_level.txt"
out_file9="$out_dir/fimo_SNAI2_SNAI2_vs_other_motif_summary_rowname_all_level.txt.txt"
out_file10="$out_dir/fimo_SNAI2_SNAI2_vs_other_motif_summary_rowname_collapse_cell.txt"
out_file11="$out_dir/fimo_SNAI2_SNAI2_vs_other_motif_summary_rowname_collapse_cell_percent.txt"
#original three columns: SNAI2 peaks that has SNAI2 motifs; SNAI2 peaks that has anohter motif; SNAI2 peaks that have two motifs; 
#use awk to add two more columns: (SNAI2 peaks that have two motifs)/(SNAI2 peaks that has SNAI2 motifs); (SNAI2 peaks that have two motifs)/(SNAI2 peaks that has another motifs)
cut -f 1 -d ' ' $out_file1 | xargs -n 3 | sed 's/ /\t/g' | awk 'BEGIN{FS=OFS="\t"}{print($0, $(NF+1)=$3/$1, $(NF+2)=$3/$2)}' >  $out_file8
paste $out_file2 $out_file8 > $out_file9
#extract the 3rd column representing SNAI2 peaks that have two motifs
paste $out_file3 <(cut -f 3 $out_file8 | xargs -n 4 | sed 's/ /\t/g') > $out_file10
#extract the 4th column representing (SNAI2 peaks that have two motifs)/(SNAI2 peaks that has SNAI2 motifs)
paste $out_file3 <(cut -f 4 $out_file8 | xargs -n 4 | sed 's/ /\t/g') > $out_file11


