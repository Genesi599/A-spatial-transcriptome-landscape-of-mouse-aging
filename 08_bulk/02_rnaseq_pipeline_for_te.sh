#!/usr/bin/bash

out=

#######################################################################################

STAR=STAR
tc=TEcount

#######################################################################################

index=
gtf=
gtf_te=

#######################################################################################

thread=

for raw_data_dir in $out/* ; do

    sample=$(basename $raw_data_dir)
    echo $sample
    sample_dir=${out}/${sample}

    qc_dir=${sample_dir}/qc_fastp

    mapping_dir=${sample_dir}/mapping_STAR_TE
    [[ -e ${mapping_dir} ]] || mkdir ${mapping_dir}

    count_dir=${sample_dir}/count_TEcount
    [[ -e ${count_dir} ]] || mkdir ${count_dir}

    pbs=${sample_dir}/${sample}.te.pbs

    bam=${mapping_dir}/${sample}_Aligned.sortedByCoord.out.bam

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

    echo "mapping \`date\` !"

    ${STAR} \\
        --runThreadN ${thread} \\
        --genomeDir ${index} \\
        --readFilesCommand zcat \\
        --readFilesIn \\
            $(ls ${qc_dir}/*_R1_trimmed.fq.gz) \\
            $(ls ${qc_dir}/*_R2_trimmed.fq.gz) \\
        --outSAMtype BAM SortedByCoordinate \\
        --outFilterType BySJout \\
        --winAnchorMultimapNmax 100 \\
        --outFilterMultimapNmax 100 \\
        --outFilterMismatchNoverLmax 0.04 \\
        --outFileNamePrefix ${mapping_dir}/${sample}_ \\
        2>&1 |
    tee ${mapping_dir}/${sample}.log

    echo "count \`date\` !"

    ${tc} \\
        -b ${bam} \\
        --sortByPos \\
        --GTF ${gtf} \\
        --TE ${gtf_te} \\
        --project ${sample} \\
        --outdir ${count_dir} \\
    2>&1 |
    tee ${count_dir}/${sample}.log

    echo "end \`date\` !"

} 2>&1 |
tee ${sample_dir}/${sample}_te_\`date '+%Y%m%d_%H%M%S'\`.log


EOF

    echo "write ${pbs} !"
    echo `date`
    qsub ${pbs}
    echo ""

done 2>&1 | 
tee ${out}/`date '+%Y%m%d_%H%M%S'`.log
