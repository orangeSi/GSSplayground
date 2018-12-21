set -vex
perl ../prepare.data.pl --list list7 --prefix out7 --outdir . --conf main.conf
cat s2.synteny.setting.conf s1.synteny.setting.conf s4.synteny.setting.conf s3.synteny.setting.conf >feature.color.label.conf7
cat  s2.to.s4.synteny.crosslink  s1.to.s2.synteny.crosslink >feature.crossing.link.conf7
perl ../plot.genome.featureCluster.pl --list list7.new --prefix out7 --outdir . --conf main.conf
