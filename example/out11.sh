set -vex
cp feature.color.label.conf10 feature.color.label.conf11

perl ../prepare.data.pl --list list11 --prefix out9 --outdir . --conf main.11.conf

cat chr14.read_mapping.setting.conf  feature.color.label.conf11 >  feature.color.label.conf11.new
cat chr14.read_mapping.crosslink >feature.crossing.link.conf11

perl ../plot.genome.featureCluster.pl --list list11.new --prefix out11 --outdir . --conf main.11.conf
