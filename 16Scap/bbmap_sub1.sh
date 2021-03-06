#!/bin/bash
#$ -S /bin/bash
#$ -N bbmap_PE300-PE
#$ -M jincheng.wang@nyumc.org
#$ -m abe
#$ -pe threaded 8
#$ -l mem_free=4G
#$ -cwd

module load java/1.8

for f in 1-QF-PE300/*_1P.fastq.gz; do sname=`basename $f .fastq.gz`; sname=${sname%_*}; time /ifs/home/wangj50/local/bbmap/bbmap.sh in=$f in2=${f/_1P/_2P} maxindel=20 minid=0.76 rcs=t killbadpairs=t pairlen=400 rescuedist=600 apd=100 outu=2-Map-PE300/${sname}_unmapped.sam outm=2-Map-PE300/${sname}_mapped.sam printunmappedcount=t covstats=2-Map-PE300/${sname}.cov.txt threads=8 -Xmx24g 2>&1 | tee 2-Map-PE300/${sname}.bbmap.log; done
