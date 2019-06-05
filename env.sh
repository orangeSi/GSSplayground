base=$(dirname `readlink -f $BASH_SOURCE`)
#echo base is $base
## for perl library
export PERL5LIB=$base/src:$base/src/Imager-1.011/lib64/perl5/:$PERL5LIB ## for perl library

# for samtool
export PATH=/zfssz5/BC_PUB/Software/03.Soft_ALL/samtools-1.7:$PATH

perl -e 'use Imager::Font'
which samtools >/dev/null

