for f in 1-QF/*_1U.fastq.gz; do sname=`basename $f .fastq.gz`; sname=${sname%_*}; time ~/local/bbmap/bbwrap.sh in=$f,${f/_1U/_2U} maxindel=20 minid=0.76 outm=3-BBMAP/${sname}_mapped_1U.sam,3-BBMAP/${sname}_mapped_2U.sam printunmappedcount=t 2>&1 | tee 3-BBMAP/${sname}_SE.bbmap.log; done

for f in 1-QF2/*_1P.fastq.gz; do sname=`basename $f .fastq.gz`; sname=${sname%_*}; time ~/local/bbmap/bbmap.sh in=$f in2=${f/_1P/_2P} maxindel=20 minid=0.76 rcs=t killbadpairs=t pairlen=400 rescuedist=600 apd=100 outm=3-BBMAP/${sname}_mapped.sam printunmappedcount=t covstats=3-BBMAP/${sname}.cov.txt 2>&1 | tee 3-BBMAP/${sname}.bbmap.log; done

for f in 1-QF2/*_1U.fastq.gz; do sname=`basename $f .fastq.gz`; sname=${sname%_*}; time ~/local/bbmap/bbwrap.sh in=$f,${f/_1U/_2U} maxindel=20 minid=0.76 outm=3-BBMAP/${sname}_mapped_1U.sam,3-BBMAP/${sname}_mapped_2U.sam printunmappedcount=t 2>&1 | tee 3-BBMAP/${sname}_SE.bbmap.log; done
