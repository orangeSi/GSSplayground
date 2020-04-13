#!/bin/sh
base=`dirname $(readlink -f $0)`/../
base=`readlink -f $base`
if [ ! -f "$base/bin/env.sh" ];
then
	echo -e "\n$(tput setaf 1)error: cannot find $base/bin/env.sh$(tput setaf 7)\n"
	exit 1
fi

#. $base/bin/env.sh $base
# check env
home=$base
#export PERL5LIB=$home/src:$home/src/Imager-1.011/lib64/perl5/:$PERL5LIB ## for perl library
export PERL5LIB=$home/src:$PERL5LIB ## for perl library

# for samtool
export PATH=/zfssz5/BC_PUB/Software/03.Soft_ALL/samtools-1.7:$PATH

perl -e 'use Imager::Font'
if [ $? -eq 0 ];
then	
	echo -e "\nfind Imager::Font"
else
	echo -e "\n$(tput setaf 1)error, not find Imager::Font, maybe you should try to re-install perl package Imager::Font(https://cpan.metacpan.org/authors/id/A/AD/ADDI/Imager-0.41.tar.gz) by cpan or manually-install$(tput setaf 7)\n"
fi
dep="samtools sort perl"
for i in $dep
do
	$i --help >/dev/null
	if [ "$?" != "0" ];
	then
		echo -e "error: not find $i"
		exit 1
	else
		j=`which $i`
		echo find $j
	fi
done




if [ $# -ne 4 ];
then
	echo -e "\nusage:\n sh $0 $(tput setaf 3)track.list$(tput setaf 7) prefix outdir/ $(tput setaf 3)main.conf$(tput setaf 7)\n\nany question, go to https://github.com/orangeSi/ClustersPloter/issues"
	exit 1
fi
#set -vex

#samtools view -h S20K.sort.bam s160:155612-178608 |samtools view - -b >S20K.s160.155612.178608.sort.bam #samtools index S20K.s160.155612.178608.sort.bam

#base=`dirname $(readlink -f $0)`/src/
#source $base/../env.sh
list=$1
prefix=$2
outdir=$3
conf=$4

for i in $list $conf
do
	if [ ! -f "$i" ];
	then
		echo -e "\n$(tput setaf 1)error: file $i not exists$(tput setaf 7)"
		exit 1
	fi
done

if [ ! -d "$outdir" ];
then
	mkdir -p $outdir
fi

cd $outdir

num=`cat $conf|grep -iE "^\s*hist_scatter_line\s*|^\s*reads_mapping\s*|^\s*synteny\s*"|wc -l`
error_flag=0
if [ "$num" -ge 1 ];
then
	echo 
	date
	echo -e "$(tput setaf 3)perl $base/src/prepare.data.pl --list $list --prefix $prefix --outdir . --conf $conf >$prefix.prepare.data.log 2>$prefix.prepare.data.error.tmp $(tput setaf 7)"
	#echo "perl $base/src/prepare.data.pl --list $list --prefix $prefix --outdir . --conf $conf >$prefix.prepare.data.log 2>$prefix.prepare.data.error.tmp"
	perl $base/src/prepare.data.pl --list $list --prefix $prefix --outdir . --conf $conf >$prefix.prepare.data.log 2>$prefix.prepare.data.error.tmp
	error_flag=$?
	#cat $prefix.prepare.data.error.tmp|grep -v '^+ ' > $prefix.prepare.data.error && rm $prefix.prepare.data.error.tmp
	cat $prefix.prepare.data.error.tmp > $prefix.prepare.data.error && rm $prefix.prepare.data.error.tmp
	echo
else
	cp $list $list.$prefix
	cp $conf $conf.$prefix
fi

if [ $error_flag -ge 1 ];
then
	error=`cat $prefix.prepare.data.error`
	echo -e "$(tput setaf 2)$error$(tput setaf 7)\n\n$(tput setaf 1)error in file: $prefix.prepare.data.error$(tput setaf 7)\n\n"
	date
	exit 1
else
date
echo
echo 
date

error_flag=0
echo -e "$(tput setaf 3)perl $base/src/plot.genome.featureCluster.pl --list $list.$prefix --prefix $prefix --outdir . --conf $conf.$prefix >$prefix.plot.log 2>$prefix.plot.error.tmp$(tput setaf 7)"
#echo "perl $base/src/plot.genome.featureCluster.pl --list $list.$prefix --prefix $prefix --outdir . --conf $conf.$prefix >$prefix.plot.log 2>$prefix.plot.error.tmp"
perl $base/src/plot.genome.featureCluster.pl --list $list.$prefix --prefix $prefix --outdir . --conf $conf.$prefix >$prefix.plot.log 2>$prefix.plot.error.tmp
error_flag=$?
cat $prefix.plot.error.tmp >$prefix.plot.error && rm $prefix.plot.error.tmp
echo
fi
echo 
#ls -ltrh $prefix.prepare.data.log $prefix.plot.log 2>/dev/null
#ls -ltrh $prefix.prepare.data.error $prefix.plot.error 2>/dev/null
if [  $error_flag -ge 1 ];
then
	error=`cat $prefix.plot.error`
	echo -e "$(tput setaf 2)$error$(tput setaf 7)\n\n$(tput setaf 1)detail error in file $prefix.plot.error, $(tput setaf 7)\n\n"
	date
	exit 1
else
	echo -e "\n\n$(tput setaf 2)finished, no error~$(tput setaf 7)\n\n"
fi
date

echo "$(tput setaf 2)output is:$(tput setaf 7)"
ls -thl $prefix.svg $prefix.notitle.svg $prefix.html
echo 
