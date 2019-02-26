#!/bin/sh
if [ $# -ne 4 ];
then
	echo -e "usage:\n sh $0 track.list prefix outdir main.conf\n\nany question, go to https://github.com/orangeSi/ClustersPloter/issues"
	exit
fi
set -vex

dep="samtools sort perl"
for i in $dep
do
	which $i 
done

#samtools view -h S20K.sort.bam s160:155612-178608 |samtools view - -b >S20K.s160.155612.178608.sort.bam #samtools index S20K.s160.155612.178608.sort.bam

base=`dirname $(readlink -f $0)`/src/
list=$1
prefix=$2
outdir=$3
conf=$4

if [ -d "$outdir" ];
then
	mkdir -p $outdir
fi

cd $outdir

num=`cat $conf|grep -E "^\s*hist_scatter_line\s*=|^\s*reads_mapping\s*=|^\s*synteny\s*="|wc -l`
if [ "$num" -ge 1 ];
then
	perl $base/prepare.data.pl --list $list --prefix $prefix --outdir . --conf $conf
else
	cp $list $list.$prefix
	cp $conf $conf.$prefix
fi


perl $base/plot.genome.featureCluster.pl --list $list.$prefix --prefix $prefix --outdir . --conf $conf.$prefix

echo $0 finish~
