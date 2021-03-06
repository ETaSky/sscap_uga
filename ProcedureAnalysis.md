#Procedure for analying 16Scap project
Created in 20171214

## For amplicon sequences
### Remove primer sequences
```shell
for f in ../../02-Data/AmpliconSeq/BG*R1.fastq.gz; do sname=${f##*/}; sname=${sname%_*}; cutadapt -j 1 -g GTGCCAGCMGCCGCGGTAA -G GACTACHVGGGTATCTAATCC -e 0.1 -n 1 -O 15 -o ${sname}_trim_R1.fastq.gz --untrimmed-output=${sname}_untrim_R1.fastq.gz -p ${sname}_trim_R2.fastq.gz --untrimmed-paired-output=${sname}_untrim_R2.fastq.gz $f ${f%R1.fastq.gz}R2.fastq.gz; done 2>&1 | tee BG_cutadapt.log

for f in ../../02-Data/AmpliconSeq/*JW-R1.fastq.gz; do sname=${f##*/}; sname=${sname%%-*}; cutadapt -j 1 -g CCTACGGGNGGCWGCAG -G GACTACHVGGGTATCTAATCC -e 0.1 -n 1 -O 15 -o JW${sname}_trim_R1.fastq.gz --untrimmed-output=JW${sname}_untrim_R1.fastq.gz -p JW${sname}_trim_R2.fastq.gz --untrimmed-paired-output=JW${sname}_untrim_R2.fastq.gz $f ${f%R1.fastq.gz}R2.fastq.gz; done 2>&1 | tee JW_cutadapt.log

cutadapt -j 1 -g CCTACGGGNGGCWGCAG -G GACTACHVGGGTATCTAATCC -e 0.1 -n 1 -O 15 -o Zymo-Amplicon_trim_R1.fastq.gz --untrimmed-output=Zymo-Amplicon_untrim_R1.fastq.gz -p Zymo-Amplicon_trim_R2.fastq.gz --untrimmed-paired-output=Zymo-Amplicon_untrim_R2.fastq.gz ../../../Tmp-Data/Zymo-Amplicon-R1.fastq.gz ../../../Tmp-Data/Zymo-Amplicon-R2.fastq.gz 2>&1 | tee Zymo_cutadapt.log 
```

### Enriched sequences
#### Quality trimming
```shell
for f in ../../02-Data/SCap/*_R1_001.fastq; do sname=`basename $f`; sname=${sname#16Scap_}; sname=${sname%%_S*}; java -jar ~/local/Trimmomatic-0.36/trimmomatic-0.36.jar PE -threads 4 -phred33 -trimlog ${sname}.log $f ${f%_R1*}_R2_001.fastq -baseout ${sname}.fastq.gz SLIDINGWINDOW:3:20 MINLEN:150; done
for f in ../../02-Data/SCap/PE300/*_R1_001.fastq.gz; do sname=`basename $f .fastq.gz`; sname=${sname#16Scap_}; sname=${sname%%_*}; java -jar ~/local/Trimmomatic-0.36/trimmomatic-0.36.jar PE -threads 4 -phred33 $f ${f/_R1/_R2} -baseout 1-QF-PE300/JT-${sname}.fastq.gz SLIDINGWINDOW:3:20 MINLEN:150; echo; done 
for f in ../../02-Data/SCap/PE150/*_R1*.fastq.gz; do sname=`basename $f .fastq.gz`; sname=${sname#16Scap_}; sname=${sname%%_*}; java -jar ~/local/Trimmomatic-0.36/trimmomatic-0.36.jar PE -threads 4 -phred33 $f ${f/_R1/_R2} -baseout 1-QF-PE150/${sname}.fastq.gz SLIDINGWINDOW:3:20 MINLEN:75; echo; done 
for f in ../../02-Data/SCap/PE150/*_R1*.fastq.gz; do sname=`basename $f .fastq.gz`; sname=${sname#16Scap_}; sname=${sname%%_*}-mock; java -jar ~/local/Trimmomatic-0.36/trimmomatic-0.36.jar PE -threads 4 -phred33 $f ${f/_R1/_R2} -baseout 1-QF-PE150/${sname}.fastq.gz SLIDINGWINDOW:3:20 MINLEN:75; echo; done 
```
#### Mapping to Greengene 13_8
```shell
## Prepare the index for the database, copy the 97_otus.fasta file from QIIME out. You can know where the file is by print_qiime_config.py in QIIME.
bbmap.sh ref=97_otus.fasta
## PE150
### Align both pair passed QF
for f in 1-QF-PE150/*_1P.fastq.gz; do sname=`basename $f .fastq.gz`; sname=${sname%_*}; time ~/local/bbmap/bbmap.sh in=$f in2=${f/_1P/_2P} maxindel=20 minid=0.76 rcs=t killbadpairs=t pairlen=400 rescuedist=600 apd=100 outm=2-Map-PE150/${sname}_mapped.sam printunmappedcount=t covstats=2-Map-PE150/${sname}.cov.txt 2>&1 | tee 2-Map-PE150/${sname}.bbmap.log; done
### Align single end
for f in 1-QF-PE150/*_1U.fastq.gz; do sname=`basename $f .fastq.gz`; sname=${sname%_*}; time ~/local/bbmap/bbwrap.sh in=$f,${f/_1U/_2U} maxindel=20 minid=0.76 outm=2-Map-PE150/${sname}_mapped_1U.sam,2-Map-PE150/${sname}_mapped_2U.sam printunmappedcount=t 2>&1 | tee 2-Map-PE150/${sname}_SE.bbmap.log; done
## PE300
### Align both pair passed QF
for f in 1-QF-PE300/*_1P.fastq.gz; do sname=`basename $f .fastq.gz`; sname=${sname%_*}; time /ifs/home/wangj50/local/bbmap/bbmap.sh in=$f in2=${f/_1P/_2P} maxindel=20 minid=0.76 rcs=t killbadpairs=t pairlen=400 rescuedist=600 apd=100 outu=2-Map-PE300/${sname}_unmapped.sam outm=2-Map-PE300/${sname}_mapped.sam printunmappedcount=t covstats=2-Map-PE300/${sname}.cov.txt threads=8 -Xmx24g 2>&1 | tee 2-Map-PE300/${sname}.bbmap.log; done
### Align single end
for f in 1-QF-PE300/*_1U.fastq.gz; do sname=`basename $f .fastq.gz`; sname=${sname%_*}; time ~/local/bbmap/bbwrap.sh in=$f,${f/_1U/_2U} maxindel=20 minid=0.76 outm=2-Map-PE300/${sname}_mapped_1U.sam,2-Map-PE300/${sname}_mapped_2U.sam printunmappedcount=t threads=auto 2>&1 | tee 2-Map-PE300/${sname}_SE.bbmap.log; done
```
#### Construct OTU table
For reads that both ends passed QF, differentiate those both ends mapped to the same ref, only one end mapped to a ref, or two ends mapped to different ref
Both ends map to the same ref, only one end was needed, so use SAM flag "83" and "99":
`for f in 2-Map-PE150/*_mapped.sam; do grep -v "^@" $f | awk 'BEGIN{FS="\t";OFS="\t"}{if($2=="83" || $2=="99") print $0}' > ${f%_mapped.sam}_both.txt; mv ${f%_mapped.sam}_both.sam 3-Filter-PE150/; echo "Finished $f"; done`
`for f in 2-Map-PE300/*_mapped.sam; do grep -v "^@" $f | awk 'BEGIN{FS="\t";OFS="\t"}{if($2=="83" || $2=="99") print $0}' > ${f%_mapped.sam}_both.sam; mv ${f%_mapped.sam}_both.sam 3-Filter-PE300/ ;echo "Finished $f"; done`
One end to a ref, only the end that mapped was needed, so use SAM flag "73","89","137", and "153":
```shell
for f in 2-Map-PE150/*_mapped_1U.sam; do grep -v "^@" $f | awk 'BEGIN{FS="\t";OFS="\t"}{if($2=="89" || $2=="73" || $2=="137" || $2=="153") print $0}' > ${f%.sam}_filter.sam ; mv ${f%.sam}_filter.sam 3-Filter-PE150/ ; done
for f in 2-Map-PE150/*_mapped_2U.sam; do grep -v "^@" $f | awk 'BEGIN{FS="\t";OFS="\t"}{if($2=="89" || $2=="73" || $2=="137" || $2=="153") print $0}' > ${f%.sam}_filter.sam ; mv ${f%.sam}_filter.sam 3-Filter-PE150/ ; done
for f in 2-Map-PE300/*_mapped_1U.sam; do grep -v "^@" $f | awk 'BEGIN{FS="\t";OFS="\t"}{if($2=="89" || $2=="73" || $2=="137" || $2=="153") print $0}' > ${f%.sam}_filter.sam ; mv ${f%.sam}_filter.sam 3-Filter-PE300/ ; done
for f in 2-Map-PE300/*_mapped_2U.sam; do grep -v "^@" $f | awk 'BEGIN{FS="\t";OFS="\t"}{if($2=="89" || $2=="73" || $2=="137" || $2=="153") print $0}' > ${f%.sam}_filter.sam ; mv ${f%.sam}_filter.sam 3-Filter-PE300/ ; done
```
#### Calculate OTU counts
This command will produce two columns, 1st col is the ID of the OTU, 2nd col is the counts. One can use the 1st col and the file 97_otu_taxonomy.txt in QIIME to extract taxonomic information for each OTU
```shell
for f in 3-Filter-PE150/*_both.sam; do sname=`basename $f _both.sam`; cat ${f%both.sam}* | awk 'BEGIN{FS="\t";OFS="\t"}{OTU[$3]+=1}END{for (i in OTU) print i,OTU[i]}' > 4-OTU-PE150/${sname}_OTU_ct.txt; done
for f in 3-Filter-PE300/*_both.sam; do sname=`basename $f _both.sam`; cat ${f%both.sam}* | awk 'BEGIN{FS="\t";OFS="\t"}{OTU[$3]+=1}END{for (i in OTU) print i,OTU[i]}' > 4-OTU-PE300/${sname}_OTU_ct.txt; done
```
Combine the OTU counts together
```shell
for f in 4-OTU-PE150/*_ct.txt; do awk -v sname="`basename $f _OTU_ct.txt`" 'BEGIN{FS="\t";OFS="\t"}{print $0, sname, "PE150"}' $f; done > Enriched_PE150_OTU_count.txt
for f in 4-OTU-PE300/*_ct.txt; do awk -v sname="`basename $f _OTU_ct.txt`" 'BEGIN{FS="\t";OFS="\t"}{print $0, sname, "PE300"}' $f; done > Enriched_PE300_OTU_count.txt

```
### Unenriched sequences
#### Combining together
```shell
for f in 1strun/16Scap*R1*.fastq.gz; do sname=`basename $f _R1_001.fastq.gz`; cat `find . -name "${sname}_*R1*.fastq.gz"` > ${sname#16Scap_}_R1.fastq.gz; done
for f in 1strun/16Scap*R2*.fastq.gz; do sname=`basename $f _R2_001.fastq.gz`; cat `find . -name "${sname}_*R2*.fastq.gz"` > ${sname#16Scap_}_R2.fastq.gz; done
```
#### Quality trimming
Using phoenix cluster
```shell
# Job submission script

#!/bin/bash
#$ -S /bin/bash
#$ -N QF
#$ -M jincheng.wang@nyumc.org
#$ -m abe
#$ -pe threaded 12
#$ -l mem_free=4G
#$ -cwd

module load trimmomatic/0.33

for f in 0-Data/*_R1.fastq.gz; do sname=`basename $f _R1.fastq.gz`; java -jar /ifs/home/wangj50/local/Trimmomatic-0.36/trimmomatic-0.36.jar PE -phred33 -threads 12 $f ${f/_R1/_R2} -baseout 1-QF/${sname}.fastq.gz SLIDINGWINDOW:3:20 MINLEN:75; done
```

#### Metaphlan
```shell
# Job Submission script

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

for f in 1-QF/*_1P.fastq.gz
do
    sname=`basename $f _1P.fastq.gz`
    metaphlan2.py $f,${f/_1P/_2P},${f/_1P/_1U},${f/_1P/_2U} --input_type fastq  --bowtie2out 2-Metaphlan/${sname}.bowtie2.bz2 --nproc 15 -t 'rel_ab_w_read_stats'  > 2-Metaphlan/${sname}_profile.txt
done
##### Combine the reads
```
for f in *profile.txt; do awk 'BEGIN{FS="\t";OFS="\t"}{print $1,$2}' $f > ${f/.txt/-rel.txt}; done
