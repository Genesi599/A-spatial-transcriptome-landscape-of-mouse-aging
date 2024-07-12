#!/usr/bin/bash

data_dir=
out=

#######################################################################################

fastp=fastp
hisat2=hisat2
samtools=samtools
fc=featureCounts

#######################################################################################

index=
gtf=

#######################################################################################

thread=10

for raw_data_dir in $data_dir/* ; do

    sample=$(basename $raw_data_dir)
    echo $sample
    sample_dir=${out}/${sample}
    [[ -e ${sample_dir} ]] || mkdir ${sample_dir}

    qc_dir=${sample_dir}/qc_fastp
    [[ -e ${qc_dir} ]] || mkdir ${qc_dir}

    mapping_dir=${sample_dir}/mapping_hisat2
    [[ -e ${mapping_dir} ]] || mkdir ${mapping_dir}

    count_dir=${sample_dir}/count_featurecounts
    [[ -e ${count_dir} ]] || mkdir ${count_dir}

    pbs=${sample_dir}/${sample}.pbs

    bam_sort=${mapping_dir}/${sample}.hisat2.sorted.bam
    bam_sort_index=${bam_sort}.bai

cat >${pbs} <<-EOF
#!/bin/bash

#PBS -N ${sample}
#PBS -o ${sample_dir}/${sample}.out
#PBS -e ${sample_dir}/${sample}.err
#PBS -d ${sample_dir}
#PBS -l nodes=1:ppn=${thread}
#PBS -q long
#PBS -V

{

    echo "qc \`date\` !"

    ${fastp} \\
        -w ${thread} \\
        -i $raw_data_dir/${sample}_R1.fq.gz \\
        -I $raw_data_dir/${sample}_R2.fq.gz \\
        -o ${qc_dir}/${sample}_R1_trimmed.fq.gz \\
        -O ${qc_dir}/${sample}_R2_trimmed.fq.gz \\
        --report_title ${sample} \\
        --html ${qc_dir}/${sample}.html \\
        --json ${qc_dir}/${sample}.json \\
        2>&1 |
    tee ${qc_dir}/${sample}.log

    echo "mapping \`date\` !"

    ${hisat2} \\
        -p ${thread} \\
        -x ${index} \\
        -1 ${qc_dir}/${sample}_R1_trimmed.fq.gz \\
        -2 ${qc_dir}/${sample}_R2_trimmed.fq.gz \\
        --novel-splicesite-outfile ${mapping_dir}/${sample}.junctions \\
        2>${mapping_dir}/${sample}.log |
    ${samtools} sort -@ ${thread} -o ${bam_sort}
    ${samtools} index -@ ${thread} ${bam_sort} ${bam_sort_index}

    echo "count \`date\` !"

    ${fc} \\
        -T ${thread} \\
        -a ${gtf} \\
        -o ${count_dir}/${sample}.featureCounts.count \\
        -p \\
        ${bam_sort} \\
        2>&1 |
    tee ${count_dir}/${sample}.log

    echo "end \`date\` !"

} 2>&1 |
tee ${sample_dir}/${sample}_\`date '+%Y%m%d_%H%M%S'\`.log


EOF

    echo "write ${pbs} !"
    echo `date`
    qsub ${pbs}
    echo ""

done 2>&1 | 
tee ${out}/`date '+%Y%m%d_%H%M%S'`.log

