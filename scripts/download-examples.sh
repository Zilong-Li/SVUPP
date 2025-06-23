#!/usr/bin/env bash

outdir="tests"
echo "download reference fasta"
if [ ! -s $outdir/1KG_ONT_VIENNA_hg38.fa ];then
    wget -N -r --no-parent --no-directories --directory-prefix=$outdir https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1KG_ONT_VIENNA/reference/1KG_ONT_VIENNA_hg38.fa.gz

    echo "gunzip fa.gz and samtools faidx"
    gunzip $outdir/1KG_ONT_VIENNA_hg38.fa.gz
    samtools faidx $outdir/1KG_ONT_VIENNA_hg38.fa
fi

echo "download aligned bam/cram files"
wget -N -r --no-parent --no-directories --directory-prefix=$outdir https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1KG_ONT_VIENNA/hg38/NA12878.hg38.cram
wget -N -r --no-parent --no-directories --directory-prefix=$outdir https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1KG_ONT_VIENNA/hg38/NA12878.hg38.cram.crai

echo "download phased reference panel"

for chr in chr{21..22};do
    wget -N -r --no-parent --no-directories --directory-prefix=$outdir https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000G_2504_high_coverage/working/20201028_3202_phased/CCDG_14151_B01_GRM_WGS_2020-08-05_${chr}.filtered.shapeit2-duohmm-phased.vcf.gz
    wget -N -r --no-parent --no-directories --directory-prefix=$outdir https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000G_2504_high_coverage/working/20201028_3202_phased/CCDG_14151_B01_GRM_WGS_2020-08-05_${chr}.filtered.shapeit2-duohmm-phased.vcf.gz.tbi
    wget -N -r --no-parent --no-directories --directory-prefix=$outdir https://raw.githubusercontent.com/rwdavies/QUILT/refs/heads/master/maps/hg38/CEU-${chr}-final.b38.txt.gz
done

echo "download SVlist"
wget -N -r --no-parent --no-directories --directory-prefix=$outdir https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1KG_ONT_VIENNA/release/v1.1/graph-augmentation/delins.sniffles.hg38.liftedT2T.13Nov2023.nygc.vcf.gz
wget -N -r --no-parent --no-directories --directory-prefix=$outdir https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1KG_ONT_VIENNA/release/v1.1/graph-augmentation/delins.sniffles.hg38.liftedT2T.13Nov2023.nygc.vcf.gz.tbi
