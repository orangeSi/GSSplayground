#!/bin/sh
if [ $# -ne 4 ];
then
	echo -e "usage:\n sh $0 track.list prefix outdir main.conf\n\nany question, go to https://github.com/orangeSi/ClustersPloter/issues"
	exit
fi
#set -vex

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
	date
	cmd="perl $base/prepare.data.pl --list $list --prefix $prefix --outdir . --conf $conf"
	echo $cmd
	perl $base/prepare.data.pl --list $list --prefix $prefix --outdir . --conf $conf >$prefix.prepare.data.log 2>$prefix.prepare.data.error.tmp
	cat $prefix.prepare.data.error.tmp|grep -v '^+ ' > $prefix.prepare.data.error && rm $prefix.prepare.data.error.tmp
	date
else
	cp $list $list.$prefix
	cp $conf $conf.$prefix
fi

if [ -s $prefix.prepare.data.error ];
then
	echo -e "\n\nerror: $prefix.prepare.data.error\n\n"
	exit
else

date
cmd="perl $base/plot.genome.featureCluster.pl --list $list.$prefix --prefix $prefix --outdir . --conf $conf.$prefix "
echo $cmd
perl $base/plot.genome.featureCluster.pl --list $list.$prefix --prefix $prefix --outdir . --conf $conf.$prefix >$prefix.plot.log 2>$prefix.plot.error.tmp
cat $prefix.plot.error.tmp|grep -v '^+ ' >$prefix.plot.error && rm $prefix.plot.error.tmp
date

fi

ls -l $prefix.prepare.data.log $prefix.plot.log $prefix.prepare.data.error $prefix.plot.error

if [ -s $prefix.plot.error ];
then
	echo -e "\n\nerror: $prefix.plot.error, $cmd\n\n"
else
	echo -e "\n\nfinish, no error~\n\n"
	ls -tl $prefix.svg $prefix.notitle.svg $prefix.html
fi
