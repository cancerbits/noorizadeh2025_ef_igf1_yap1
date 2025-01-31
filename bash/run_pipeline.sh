#!/bin/bash

source parse_yaml.sh
eval $(parse_yaml ../config.yaml CONF_)

mkdir -p $CONF_project_out_host/submission/ 
mkdir -p $CONF_project_out_host/pipeline/
 
 
 
export PEPENV="$CONF_metadata_pipeline_interfaces"

nohup looper rerun --sp atac $CODEBASE/side_projects/$PROJECT/metadata/config.yaml > "$CODEBASE/side_projects/$PROJECT/nohup_$(hostname)_${USER}_atac.log" 2>&1 &

# run additional TSS score calculation for ATAC-seq samples (this step was skipped in the initial pipeline run because of a missing TSS annotation file):
for S in $(cut -f2 -d, $CODEBASE/side_projects/$PROJECT/metadata/samples_atac.csv) ; do
	INFILE=$OUT/ewing_mouse/results/pipeline/${S}/aligned_mm10/${S}_sort_dedup.bam
	if [ -f $INFILE ] ; then
		echo $S
		OUTFILE=$OUT/ewing_mouse/results/pipeline/${S}/QC_mm10/${S}_TssEnrichment.txt
		if [ ! -f $OUTFILE ] ; then
			/data_synology_rg3/cancerbits/tools/pepatac/tools/pyTssEnrichment.py -a $INFILE -b ${RESOURCES}/genomes/refgenie/mm10/mm10_TSS.tsv -p ends -c 8 -e 2000 -u -v -s 4 -o $OUTFILE
		fi
	fi
done

