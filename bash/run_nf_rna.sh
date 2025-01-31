#!/bin/bash

source parse_yaml.sh
eval $(parse_yaml ../config.yaml CONF_)

DATAHOME=${CONF_data_root_host}/
OUTHOME=${CONF_out_root_host}/pipeline/nf-core/rna
FQ_DIR=${OUTHOME}/fastq
SASH=${CONF_project_root_host}/metadata/samples_nf_rna.csv
SASH_NF=${OUTHOME}/sash_rnaseq.csv
SAMTOOLS=${CONF_samtools}
PIPE_NAME=nf-core/rnaseq
PIPE_VER=3.14.0
export JAVA_HOME=${CONF_java_home}
export JAVA_CMD=${JAVA_HOME}/bin/java
export PATH=${CONF_nextflow_home}:${JAVA_HOME}/bin:${PATH}
export NXF_VER=23.10.1
export NXF_SINGULARITY_CACHEDIR=${OUTHOME}/singu_cache/
mkdir -p ${NXF_SINGULARITY_CACHEDIR}

ORG=${CONF_implied_attributes_organism_mouse_genome}
GTF=${CONF_resources_root_host}/genomes/${ORG}/${ORG}.gtf
FA=${CONF_resources_root_host}/genomes/${ORG}/${ORG}.fa


mkdir -p ${FQ_DIR}

cd ${OUTHOME}



echo "Prepare input files..."

echo "sample,fastq_1,fastq_2,strandedness" > ${SASH_NF}
while IFS=',' read -r SAMPLE BAM REST ; do 
	if [ -f ${DATAHOME}/${BAM} ] ; then
		echo "  ${SAMPLE}"
		
		FQ="${FQ_DIR}/${SAMPLE}.fastq.gz"
		
		if [ ! -f ${FQ} ] ; then
			echo "	  - $BAM --> $FQ"
			${SAMTOOLS} fastq --threads ${CONF_nextflow_cpus} ${DATAHOME}/${BAM} | gzip > ${FQ}
		fi
		
		echo "${SAMPLE},${FQ},,auto" >> ${SASH_NF}
	fi
done < ${SASH}



echo "Config:"
echo "  Organism: $ORG"
echo "  GTF: $GTF"
echo "  FASTA: $FA"

echo "Run nextflow $PIPE_NAME ($PIPE_VER)..."

nextflow run ${PIPE_NAME} \
	-r ${PIPE_VER} \
	-resume \
	-profile singularity \
	--max_cpus ${CONF_nextflow_cpus} \
	--max_memory ${CONF_nextflow_mem} \
	-work-dir ${OUTHOME}/work \
    --fasta ${FA} \
    --gtf ${GTF} \
    --gencode \
    --aligner 'star_rsem' \
    --save_reference \
	--input ${SASH_NF} \
    --outdir ${OUTHOME}

echo "DONE!"

