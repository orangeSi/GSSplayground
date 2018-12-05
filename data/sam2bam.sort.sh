
bwa=/zfssz5/BC_PUB/Software/03.Soft_ALL/bwa-0.7.17/bwa
samtools=/zfssz5/BC_PUB/Software/03.Soft_ALL/samtools-1.7/samtools

if [ "$2"  == "" ];
then
	tmp=`echo $samtools|awk '{print $NF}'`
	for i in  $tmp
	do
		if [ ! -f "$i" ];
		then
			echo "check error: $i not exists"
			exit
		fi
	done
	echo "sh $0 <sam> <prefix>"
	exit
fi
sam=$1
bam=$2


set -vex
$samtools view -bS $sam > ${sam}.bam 
$samtools sort ${sam}.bam -o ${bam}.sort.bam
$samtools index ${bam}.sort.bam
rm ${sam}.bam $sam

