
set -vex
r1=s3.seq.bwa.read1.fastq.gz
r2=s3.seq.bwa.read2.fastq.gz
/zfssz5/BC_PUB/Software/03.Soft_ALL/bwa-0.7.17/bwa index s3.seq
/zfssz5/BC_PUB/Software/03.Soft_ALL/bwa-0.7.17/bwa mem s3.seq s3.seq.bwa.read1.fastq.gz -t 3 >./s3.seq.sam
sam=./s3.seq.sam
/zfssz5/BC_PUB/Software/03.Soft_ALL/samtools-1.7/samtools view -bS $sam > $sam.bam 
/zfssz5/BC_PUB/Software/03.Soft_ALL/samtools-1.7/samtools sort $sam.bam -o $sam.sorted.bam
rm $sam.bam $sam
touch ./s3.seq.bwa.sh.sign

