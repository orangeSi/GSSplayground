set -vex
perl ../plot.genome.featureCluster.pl --list list5 --prefix out --outdir . --conf main.conf

# remove the crosslinks and decrease the height of figure, change feature type to rect
perl ../plot.genome.featureCluster.pl --list list5 --prefix out2 --outdir . --conf main.2.conf

perl ../plot.genome.featureCluster.pl --list list5 --prefix out3 --outdir . --conf main.3.conf
#perl ../plot.genome.featureCluster.pl --list list5.fake --prefix out.fake --outdir . --conf main.fake.conf # not support yet
perl ../plot.genome.featureCluster.pl --list list6 --prefix out6 --outdir . --conf main.6.conf
perl ../plot.genome.featureCluster.pl --list list7 --prefix out7 --outdir . --conf main.conf
