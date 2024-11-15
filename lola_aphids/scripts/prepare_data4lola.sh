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
# agat_convert_sp_gff2bed.pl --gff ${sample}_OR.gff --sub gene -c ../../scripts/agat_config.yaml --out ../bed/${sample}_gene_OR.bed
# agat_convert_sp_gff2bed.pl --gff ${sample}_GR.gff --sub gene -c ../../scripts/agat_config.yaml --out ../bed/${sample}_gene_GR.bed

cd $path/data/bed
bedtools slop -i ${sample}_gene.bed -g ../genomes/${sample}.fasta.fai -b 2000 > ${sample}_2kb.bed
bedtools slop -i ${sample}_gene_OR.bed -g ../genomes/${sample}.fasta.fai -b 2000 > ${sample}_OR_2kb.bed
bedtools slop -i ${sample}_gene_GR.bed -g ../genomes/${sample}.fasta.fai -b 2000 > ${sample}_GR_2kb.bed

bedtools slop -i ${sample}_gene.bed -g ../genomes/${sample}.fasta.fai -b 10000 > ${sample}_10kb.bed
bedtools slop -i ${sample}_gene_OR.bed -g ../genomes/${sample}.fasta.fai -b 10000 > ${sample}_OR_10kb.bed
bedtools slop -i ${sample}_gene_GR.bed -g ../genomes/${sample}.fasta.fai -b 10000 > ${sample}_GR_10kb.bed

#perl -ane '{($symbol) = /symbol=([^;]+)/ ; print join ("\t", $F[0],$F[3], $F[4],$symbol), "\n" if $F[2] eq "gene"}' ${sample}.gff > ../bed/${sample}.bed
#perl -ane '{
#    print join ("\t", $F[0],$F[3], $F[4]), "\n" if $F[2] eq "gene";
#}' ${sample}.gff > ../bed/${sample}.bed
#perl -ane '{($symbol) = /symbol=([^;]+)/ ; print join ("\t", $F[0],$F[3], $F[4],$symbol), "\n"}' ${sample}_GR.gff > ../bed/${sample}_GR.bed
#perl -ane '{($symbol) = /symbol=([^;]+)/ ; print join ("\t", $F[0],$F[3], $F[4],$symbol), "\n"}' ${sample}_OR.gff > ../bed/${sample}_OR.bed

#perl -ane '{print join ("\t", $F[0],$F[3], $F[4]), "\n"}' ${sample}.gff3 > ${sample}.bed
#perl -ane '{print join ("\t", $F[0],$F[3], $F[4]), "\n"}' ${sample}_SFBB_only.gff3 > ${sample}_SFBB_only.bed

#perl $path/scripts/enlargebed.pl -bed ${sample}.bed -bound 2000 -genome ../genomes/${sample}.fasta.fai > ${sample}_2kb.bed
#perl $path/scripts/enlargebed.pl -bed ${sample}_GR.bed -bound 2000 -genome ../genomes/${sample}.fasta.fai > ${sample}_GR_2kb.bed
#perl $path/scripts/enlargebed.pl -bed ${sample}_OR.bed -bound 2000 -genome ../genomes/${sample}.fasta.fai > ${sample}_OR_2kb.bed

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

#for i in {1..1000}; do bedtools shuffle -i ${sample}_2000.sorted.merged.bed -g ${sample}.fa.fai >> ${sample}_2000.shuffle1000.bed ; done 
#for i in {1..1000}; do bedtools shuffle -i ${sample}_OR_2000.sorted.merged.bed -g ${sample}.fa.fai >> ${sample}_OR_2000.shuffle1000.bed ; done
#for i in {1..1000}; do bedtools shuffle -i ${sample}_OR_2000.sorted.merged.bed -g ${sample}.fa.fai >> ${sample}_OR_2000_2.shuffle1000.bed ; done
#for i in {1..1000}; do bedtools shuffle -i ${sample}_OR_2000.sorted.merged.bed -g ${sample}.fa.fai >> ${sample}_OR_2000_3.shuffle1000.bed ; done
#for i in {1..1000}; do bedtools shuffle -i ${sample}_GR_2000.sorted.merged.bed -g ${sample}.fa.fai >> ${sample}_GR_2000.shuffle1000.bed ; done
#for i in {1..1000}; do bedtools shuffle -i ${sample}_GR_2000.sorted.merged.bed -g ${sample}.fa.fai >> ${sample}_GR_2000_2.shuffle1000.bed ; done
#for i in {1..1000}; do bedtools shuffle -i ${sample}_GR_2000.sorted.merged.bed -g ${sample}.fa.fai >> ${sample}_GR_2000_3.shuffle1000.bed ; done

cd $path/data/TE_aphid_URGI/Annotation_genome_puceron/${sample}

# cat ${sample}_TEannotGr2_GFF3chr/*.gff3 > ${sample}_TE_annotation.gff3
perl $path/scripts/REPET2LOLA.pl -gff ${sample}_TE_annotation.gff3

done

