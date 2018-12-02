set -vx
perl ../prepare.data.pl --list list5 --prefix out9 --outdir . --conf main.10.conf
cat s2.setting.conf s3.setting.conf feature.color.label.conf10 >  feature.color.label.conf10.new
perl ../plot.genome.featureCluster.pl --list list5.ytick --prefix out10 --outdir . --conf main.10.conf
