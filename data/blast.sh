que="NC_017659.1.indel.fa"
ref="NC_017659.1.fa"

/ifs4/BC_PUB/biosoft/pipeline/Package/myth/bin/blast/blastn.sh $ref $que  $que.to.$ref.blast 3 100
/ifs4/BC_PUB/biosoft/pipeline/Package/myth/bin/blast/blastn.sh $ref $ref  $ref.to.$ref.blast 3 100
