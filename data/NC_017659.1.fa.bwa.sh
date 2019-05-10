
set -vex
r1=NC_017659.1.fa.bwa.read1.fastq.gz
r2=NC_017659.1.fa.bwa.read2.fastq.gz
/zfssz5/BC_PUB/Software/03.Soft_ALL/bwa-0.7.17/bwa index NC_017659.1.fa
/zfssz5/BC_PUB/Software/03.Soft_ALL/bwa-0.7.17/bwa mem NC_017659.1.fa NC_017659.1.fa.bwa.read1.fastq.gz NC_017659.1.fa.bwa.read2.fastq.gz -t 3 >./NC_017659.1.fa.sam
sam=./NC_017659.1.fa.sam
/zfssz5/BC_PUB/Software/03.Soft_ALL/samtools-1.7/samtools view -bS $sam > $sam.bam 
/zfssz5/BC_PUB/Software/03.Soft_ALL/samtools-1.7/samtools sort $sam.bam -o $sam.sorted.bam
/zfssz5/BC_PUB/Software/03.Soft_ALL/samtools-1.7/samtools index $sam.sorted.bam
rm $sam.bam $sam
touch ./NC_017659.1.fa.bwa.sh.sign

