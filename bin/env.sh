home=$(dirname `readlink -f $BASH_SOURCE`)
home="$home/../"
#echo base is $base
## for perl library
export PERL5LIB=$home/src:$home/src/Imager-1.011/lib64/perl5/:$PERL5LIB ## for perl library

# for samtool
export PATH=/zfssz5/BC_PUB/Software/03.Soft_ALL/samtools-1.7:$PATH


function check_env(){
	perl -e 'use Imager::Font'
	if [ "$?" == "0" ];
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
			exit
		else
			j=`which $i`
			echo find $j
		fi
	done
}

check_env
