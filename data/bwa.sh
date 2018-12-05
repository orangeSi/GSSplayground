
bwa=/zfssz5/BC_PUB/Software/03.Soft_ALL/bwa-0.7.17/bwa
samtools=/zfssz5/BC_PUB/Software/03.Soft_ALL/samtools-1.7/samtools

if [ "$5"  == "" ];
then
	tmp=`echo $samtools|awk '{print $NF}'`
	for i in  $bwa $tmp
	do
		if [ ! -f "$i" ];
		then
			echo "check error: $i not exists"
			exit
		fi
	done
	echo "sh $0 <r1,r2> <ref.fa>  <outdir> <prefix> <bwa cpu>"
	exit
fi
r1r2=$1
ass=$2
outdir=$3
prefix=$4
cpu=$5

r1=`echo $r1r2|awk -F ',' '{print $1}'`
r2=`echo $r1r2|awk -F ',' '{print $2}'`

if [ ! -d "$outdir" ];
then
	mkdir -p $outdir
fi
ref=$ass
out=$outdir/$prefix

echo "
set -vex
r1=$r1
r2=$r2
$bwa index $ref
$bwa mem $ref $r1 -t $cpu >$out.sam
sam=$out.sam
$samtools view -bS \$sam > \$sam.bam 
$samtools sort \$sam.bam -o \$sam.sorted.bam
rm \$sam.bam \$sam
touch $outdir/${prefix}.bwa.sh.sign
" >$outdir/${prefix}.bwa.sh

echo $outdir/${prefix}.bwa.sh wait to run 
