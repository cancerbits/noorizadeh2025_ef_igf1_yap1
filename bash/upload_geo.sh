#!/bin/bash

source parse_yaml.sh
eval $(parse_yaml ../config.yaml CONF_)

shopt -s nullglob



# RNA-seq

echo "RNA-seq:"

mkdir -p ${CONF_out_root}/geo_upload/noorizadeh_rna
while IFS="," read -r sample_name bam_file x; do 
  in=${CONF_data_root}/${bam_file}
  out=${CONF_out_root}/geo_upload/noorizadeh_rna/S_${sample_name}.bam
  if [ -f ${in} -a ! -f ${out} ]; then
    echo "- ${out}"
    cp ${in} ${out}
  fi
  
  for bw in "${CONF_out_root}/pipeline/nf-core/rna/star_rsem/bigwig/${sample_name}".*.bigWig ; do
    out=${CONF_out_root}/geo_upload/noorizadeh_rna/S_$(basename $bw)
    if [ ! -f ${out} ]; then
      echo "- ${out}"
      cp ${bw} ${out}
    fi
  done
done < ${CONF_project_root}/metadata/samples_nf_rna.csv
 
cat ${CONF_out_root}/pipeline/nf-core/rna/star_rsem/rsem.merged.gene_counts.tsv | gzip > ${CONF_out_root}/geo_upload/noorizadeh_rna/rsem.merged.gene_counts.tsv.gz
 
echo "- MD5 sums"
cd ${CONF_out_root}/geo_upload/noorizadeh_rna/
md5sum *.bam *.tsv.gz *.bigWig > md5.txt


# ATAC-seq

echo "ATAC-seq:"
 
mkdir -p ${CONF_out_root}/geo_upload/noorizadeh_atac
while IFS="," read -r library_pool sample_name use genotype mouse features description organism cell_type sample_group_detail sample_group mouse_old bio_rep bio_rep_num read_type bsf_folder flowcell lane x; do
  in=${CONF_data_root}/${bsf_folder}/${flowcell}_L${lane}/${flowcell}_${lane}#${sample_name}.bam
  out=${CONF_out_root}/geo_upload/noorizadeh_atac/S_${sample_name}.bam
  if [ -f ${in} -a ! -f ${out} ]; then
    echo "- ${out}"
    cp ${in} ${out}
  fi
   
  in=${CONF_out_root}/pipeline/${sample_name}/aligned_${CONF_implied_attributes_organism_mouse_refgenie_genome}/${sample_name}_smooth.bw
  out=${CONF_out_root}/geo_upload/noorizadeh_atac/S_${sample_name}_smooth.bw
  if [ -f ${in} -a ! -f ${out} ]; then
    echo "- ${out}"
    cp ${in} ${out}
  fi
   
  in=${CONF_out_root}/pipeline/${sample_name}/peak_calling_${CONF_implied_attributes_organism_mouse_refgenie_genome}/${sample_name}_peaks.narrowPeak
  out=${CONF_out_root}/geo_upload/noorizadeh_atac/S_${sample_name}_peaks.narrowPeak
  if [ -f ${in} -a ! -f ${out} ]; then
    echo "- ${out}"
    cp ${in} ${out}
  fi
done < ${CONF_project_root}/metadata/samples_atac.csv
 
echo "- MD5 sums"
cd ${CONF_out_root}/geo_upload/noorizadeh_atac/
md5sum *.bam > md5_raw.txt
md5sum *.bw *.narrowPeak > md5_proc.txt



# GEO upload

echo "GEO uploads"

USR=geoftp
PASSWD=XXXXX
FTP_HOST=ftp-private.ncbi.nlm.nih.gov
REMOTE_ROOT=XXXXX

echo "GEO upload RNA"
ncftpput -F -R -m -z -u $USR -p "$PASSWD" $FTP_HOST ./$REMOTE_ROOT/noorizadeh_rna ${CONF_out_root_host}/geo_upload/noorizadeh_rna/*

echo "GEO upload ATAC"

ncftpput -F -R -m -z -u $USR -p "$PASSWD" $FTP_HOST ./$REMOTE_ROOT/noorizadeh_atac ${CONF_out_root_host}/geo_upload/noorizadeh_atac/*



