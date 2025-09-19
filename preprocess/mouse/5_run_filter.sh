
mkdir -p Drp1_mouse_dge/merged/filter
awk 'NR==1 || $2~/Dnm1l/' Drp1_mouse_dge/merged/default/counts_transcript.txt > Drp1_mouse_dge/merged/filter/counts_transcript.txt
awk 'NR==1 || $2~/Dnm1l/' Drp1_mouse_dge/merged/default/CPM_transcript.txt > Drp1_mouse_dge/merged/filter/CPM_transcript.txt

grep 'Dnm1l' Drp1_mouse_dge/merged/default/extended_annotations.gtf > Drp1_mouse_dge/merged/filter/extended_annotations.gtf