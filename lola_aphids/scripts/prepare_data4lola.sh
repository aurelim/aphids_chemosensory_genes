#!/bin/bash
#SBATCH --mem=50G

module load bioinfo/samtools/1.19
module load bioinfo/bedtools/2.30.0
module load devel/Miniconda/Miniconda3
module load bioinfo/AGAT/1.2.0

path=/home/amesnil/work/lola_aphids

for sample in AGLY AGOS APIS DNOX DPLA DVIT ELAN MCER MPER PNIG RMAI RPAD SMIS
# for sample in APIS
do

cd $path/data/genomes/

samtools faidx ${sample}.fasta

cd $path/data/gff/

agat_convert_sp_gff2bed.pl --gff ${sample}.gff --sub Level1 --out ../bed/${sample}_gene.bed
cat ../bed/${sample}_gene.bed | grep -i "or" > ../bed/${sample}_gene_OR.bed
cat ../bed/${sample}_gene.bed | grep -i "gr" > ../bed/${sample}_gene_GR.bed

cd $path/data/bed
bedtools slop -i ${sample}_gene.bed -g ../genomes/${sample}.fasta.fai -b 2000 > ${sample}_2kb.bed
bedtools slop -i ${sample}_gene_OR.bed -g ../genomes/${sample}.fasta.fai -b 2000 > ${sample}_OR_2kb.bed
bedtools slop -i ${sample}_gene_GR.bed -g ../genomes/${sample}.fasta.fai -b 2000 > ${sample}_GR_2kb.bed

bedtools slop -i ${sample}_gene.bed -g ../genomes/${sample}.fasta.fai -b 10000 > ${sample}_10kb.bed
bedtools slop -i ${sample}_gene_OR.bed -g ../genomes/${sample}.fasta.fai -b 10000 > ${sample}_OR_10kb.bed
bedtools slop -i ${sample}_gene_GR.bed -g ../genomes/${sample}.fasta.fai -b 10000 > ${sample}_GR_10kb.bed


sort -k1,1 -k2,2n ${sample}_2kb.bed > ${sample}_2kb.sorted.bed
sort -k1,1 -k2,2n ${sample}_GR_2kb.bed > ${sample}_GR_2kb.sorted.bed
sort -k1,1 -k2,2n ${sample}_OR_2kb.bed > ${sample}_OR_2kb.sorted.bed

sort -k1,1 -k2,2n ${sample}_10kb.bed > ${sample}_10kb.sorted.bed
sort -k1,1 -k2,2n ${sample}_GR_10kb.bed > ${sample}_GR_10kb.sorted.bed
sort -k1,1 -k2,2n ${sample}_OR_10kb.bed > ${sample}_OR_10kb.sorted.bed

bedtools merge -i ${sample}_2kb.sorted.bed > ${sample}_2kb.sorted.merged.bed
bedtools merge -i ${sample}_GR_2kb.sorted.bed > ${sample}_GR_2kb.sorted.merged.bed
bedtools merge -i ${sample}_OR_2kb.sorted.bed > ${sample}_OR_2kb.sorted.merged.bed 

bedtools merge -i ${sample}_10kb.sorted.bed > ${sample}_10kb.sorted.merged.bed
bedtools merge -i ${sample}_GR_10kb.sorted.bed > ${sample}_GR_10kb.sorted.merged.bed
bedtools merge -i ${sample}_OR_10kb.sorted.bed > ${sample}_OR_10kb.sorted.merged.bed 


cd $path/data/TE_aphid_URGI/Annotation_genome_puceron/${sample}

cat ${sample}_TEannotGr2_GFF3chr/*.gff3 > ${sample}_TE_annotation.gff3
perl $path/scripts/REPET2LOLA.pl -gff ${sample}_TE_annotation.gff3

done

