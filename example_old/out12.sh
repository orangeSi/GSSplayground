set -vex
cp feature.color.label.conf12 feature.color.label.conf12.new
cp feature.crossing.link.conf9 feature.crossing.link.conf12
perl ../src/prepare.data.pl --list list5 --prefix out9 --outdir . --conf main.12.conf
cat  s2.read_mapping.setting.conf s2.hist_scatter_line.setting.conf  s3.hist_scatter_line.setting.conf feature.color.label.conf12 >  feature.color.label.conf12.new

cat s3.hist_scatter_line.crosslink feature.crossing.link.conf9.new >feature.crossing.link.conf12

perl ../src/plot.genome.featureCluster.pl --list list5.ytick --prefix out12 --outdir . --conf main.12.conf
