#perl -e 'print 30000*100/300'
#10000
set -vex
number_read=10000
for i in $(ls s*.seq)
do
	#/zfssz5/BC_PS/sikaiwei/software/read_simulater/DWGSIM/dwgsim $i $i -N $number_read -1 150 -2 150 -c 0 -S 2 -H -z 1 -o 1 -r 0.02 -R 0.20 -z 2018
	#sh bwa.sh ${i}.bwa.read1.fastq.gz,${i}.bwa.read2.fastq.gz $i . $i 3 && sh ${i}.bwa.sh
	#/zfssz4/BC_COM_P0/F18FTSECKF1389/ASPjisD/software/pbsim/pbsim-1.0.3-Linux-amd64/Linux-amd64/bin/pbsim --data-type CLR --depth 20 --length-min 1000 --length-max 10000 --seed 0 --sample-fastq /zfssz4/BC_COM_P0/F18FTSECKF1389/ASPjisD/data/pb/chloroplast/sd_0001.fasta  $i --prefix $i.test
	#cat ${i}.test*fastq >${i}.longreads.fq && rm $i.test*
	#/share/backup/sikaiwei/sikaiwei/assembly/sprai/sprai-0.9.9.19/fq2fa.pl  ${i}.longreads.fq >${i}.longreads.fa && gzip ${i}.longreads.fq
	minimap2 $i ${i}.longreads.fa -ax map-pb >${i}.longreads.map2ref.sam
	sh sam2bam.sort.sh ${i}.longreads.map2ref.sam ${i}.longreads.map2ref
	minimap2 $i ${i}.longreads.fa -cx map-pb |gzip >${i}.longreads.map2ref.paf.gz
	
done
