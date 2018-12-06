set -vex
cp feature.color.label.conf10 feature.color.label.conf10.new
cp feature.crossing.link.conf9 feature.crossing.link.conf9.new

perl ../prepare.data.pl --list list5 --prefix out9 --outdir . --conf main.10.conf

cat s[1234].{plot_depth,read_mapping}.setting.conf feature.color.label.conf10.old >  feature.color.label.conf10.new
cat s3.plot_depth.crosslink  feature.crossing.link.conf9 >feature.crossing.link.conf9.new

perl ../plot.genome.featureCluster.pl --list list5.ytick --prefix out10 --outdir . --conf main.10.conf
