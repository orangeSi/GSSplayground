set -vex
que="NC_017659.1.indel.fa"
ref="NC_017659.1.fa"
minimap2 $ref $que -cx asm5 |gzip >$que.to.$ref.paf.gz
