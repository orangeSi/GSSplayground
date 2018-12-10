set -vex
cp feature.color.label.conf10 feature.color.label.conf10.new
cp feature.crossing.link.conf9 feature.crossing.link.conf10
perl ../prepare.data.pl --list list5 --prefix out9 --outdir . --conf main.10.conf
cat  s2.read_mapping.setting.conf s2.plot_depth.setting.conf  s3.plot_depth.setting.conf feature.color.label.conf10 >  feature.color.label.conf10.new

cat s3.plot_depth.crosslink feature.crossing.link.conf9.new >feature.crossing.link.conf10

perl ../plot.genome.featureCluster.pl --list list5.ytick --prefix out10 --outdir . --conf main.10.conf
