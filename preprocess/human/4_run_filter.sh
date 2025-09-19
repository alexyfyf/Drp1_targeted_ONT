mkdir -p Drp1_human_dge/filter
awk 'NR==1 || $2~/DNM1L/' Drp1_human_dge/default/counts_transcript.txt > Drp1_human_dge/filter/counts_transcript.txt
awk 'NR==1 || $2~/DNM1L/' Drp1_human_dge/default/CPM_transcript.txt > Drp1_human_dge/filter/CPM_transcript.txt

grep 'DNM1L' Drp1_human_dge/default/extended_annotations.gtf > Drp1_human_dge/filter/extended_annotations.gtf
