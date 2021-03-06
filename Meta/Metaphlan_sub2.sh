#!/bin/bash
#$ -S /bin/bash
#$ -N Metaphlan2
#$ -M jincheng.wang@nyumc.org
#$ -m abe
#$ -pe threaded 15
#$ -l mem_free=4G
#$ -cwd

# The script was submitting in analysis base folder
module load samtools/1.3 metaphlan/2.0

for f in 1-QF2/*_1P.fastq.gz
do
    sname=`basename $f _1P.fastq.gz`
    metaphlan2.py $f,${f/_1P/_2P},${f/_1P/_1U},${f/_1P/_2U} --input_type fastq  --bowtie2out 2-Metaphlan/${sname}.bowtie2.bz2 --nproc 15 -t 'rel_ab_w_read_stats'  > 2-Metaphlan/${sname}_profile.txt
done
