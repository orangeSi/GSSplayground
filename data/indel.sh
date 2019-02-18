set -vex

ref=NC_017659.1.fa
que=NC_017659.1.indel.fa
que_long=NC_017659.1.indel.fa.longreads.fa
que_short_r1=NC_017659.1.indel.fa.bwa.read1.fastq.gz
que_short_r2=NC_017659.1.indel.fa.bwa.read2.fastq.gz

#sh bwa.sh $que_short_r1,$que_short_r2 $ref . $que.to.$ref.sr 3 && sh $que.to.$ref.sr.bwa.sh
minimap2 $ref $que_long -ax map-pb >$que.to.$ref.longreads.sam
sh sam2bam.sort.sh $que.to.$ref.longreads.sam $que.to.$ref.longreads
