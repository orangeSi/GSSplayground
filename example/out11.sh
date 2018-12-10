set -vex

echo "" > feature.color.label.conf11.new
echo "" >feature.color.label.conf11

perl ../prepare.data.pl --list list6 --prefix out9 --outdir . --conf main.11.conf

cat chr14.read_mapping.setting.conf  feature.color.label.conf11 >  feature.color.label.conf11.new

perl ../plot.genome.featureCluster.pl --list list6.new --prefix out11 --outdir . --conf main.11.conf
