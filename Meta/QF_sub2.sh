#!/bin/bash
#$ -S /bin/bash
#$ -N QF
#$ -M jincheng.wang@nyumc.org
#$ -m abe
#$ -pe threaded 12
#$ -l mem_free=4G
#$ -cwd

module load trimmomatic/0.33

for f in 0-Data/*.R1.fastq.gz; do sname=`basename $f .R1.fastq.gz`; java -jar /ifs/home/wangj50/local/Trimmomatic-0.36/trimmomatic-0.36.jar PE -phred33 -threads 12 $f ${f/.R1/.R2} -baseout 1-QF2/${sname}.fastq.gz SLIDINGWINDOW:3:20 MINLEN:75; done
