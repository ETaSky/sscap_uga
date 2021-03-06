#!/bin/bash
#$ -S /bin/bash
#$ -N bbmap_PE300-SE
#$ -M jincheng.wang@nyumc.org
#$ -m abe
#$ -pe threaded 8
#$ -l mem_free=4G
#$ -cwd

module load java/1.8

for f in 1-QF-PE300/*_1U.fastq.gz; do sname=`basename $f .fastq.gz`; sname=${sname%_*}; time /ifs/home/wangj50/local/bbmap/bbwrap.sh in=$f,${f/_1U/_2U} maxindel=20 minid=0.76 outu=2-Map-PE300/${sname}_mapped_1U.sam,2-Map-PE300/${sname}_mapped_2U.sam printunmappedcount=t threads=8 -Xmx24g 2>&1 | tee 2-Map-PE300/${sname}_SE.bbmap.log; done
