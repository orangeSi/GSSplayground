ref=NC_017659.1.fa
q=NC_017659.1.indel.fa
prefix=mummer
nucmer --mum -c 100 -p $prefix $ref $q
delta-filter $prefix.delta >$prefix.delta.filter
show-coords -rTl $prefix.delta.filter > $prefix.delta.filter.coords

