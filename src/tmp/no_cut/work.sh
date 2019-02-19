set -vex

perl ../../prepare.data.pl --list list5 --prefix out5 --outdir . --conf main.11.conf

cat aa.read_mapping.setting.conf > feature.color.label.conf
cat aa.read_mapping.crosslink    > feature.crossing.link

perl ../../plot.genome.featureCluster.pl --list list5.new --prefix out5 --outdir . --conf main.11.new.conf
