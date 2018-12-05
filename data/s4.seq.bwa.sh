
set -vex
r1=s4.seq.bwa.read1.fastq.gz
r2=s4.seq.bwa.read2.fastq.gz
/zfssz5/BC_PUB/Software/03.Soft_ALL/bwa-0.7.17/bwa index s4.seq
/zfssz5/BC_PUB/Software/03.Soft_ALL/bwa-0.7.17/bwa mem s4.seq s4.seq.bwa.read1.fastq.gz s4.seq.bwa.read2.fastq.gz -t 3 >./s4.seq.sam
sam=./s4.seq.sam
/zfssz5/BC_PUB/Software/03.Soft_ALL/samtools-1.7/samtools view -bS $sam > $sam.bam 
/zfssz5/BC_PUB/Software/03.Soft_ALL/samtools-1.7/samtools sort $sam.bam -o $sam.sorted.bam
/zfssz5/BC_PUB/Software/03.Soft_ALL/samtools-1.7/samtools index $sam.sorted.bam
rm $sam.bam $sam
touch ./s4.seq.bwa.sh.sign

