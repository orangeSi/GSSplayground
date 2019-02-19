set -vex
samtools view -h /zfssz2/BC_COM_P5/USER/lixiuli/greenalgae/00.assembly/mate/S20K.sort.bam s160:155612-178608 |samtools view - -b >S20K.s160.155612.178608.sort.bam
samtools index S20K.s160.155612.178608.sort.bam

perl ../prepare.data.pl --list list5 --prefix out5 --outdir . --conf main.11.conf

cat aa.read_mapping.setting.conf > feature.color.label.conf
cat aa.read_mapping.crosslink    > feature.crossing.link

perl ../plot.genome.featureCluster.pl --list list5.new --prefix out5 --outdir . --conf main.11.new.conf
