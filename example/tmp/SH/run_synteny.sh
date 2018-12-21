cd /hwfssz4/BC_COM_P0/F18FTSECKF1389/ASPjisD/ClustersPloter/example/tmp/synteny
/ifs4/BC_PUB/biosoft/pipeline/Package/MUMmer-3.22/nucmer /hwfssz4/BC_COM_P0/F18FTSECKF1389/ASPjisD/ClustersPloter/example/tmp/filter/Target.200.fa /hwfssz4/BC_COM_P0/F18FTSECKF1389/ASPjisD/ClustersPloter/example/tmp/filter/Query.200.fa -prefix a-b 2>/dev/null
/ifs4/BC_PUB/biosoft/pipeline/Package/MUMmer-3.22/delta-filter -1 a-b.delta > a-b.filter
/ifs4/BC_PUB/biosoft/pipeline/Package/MUMmer-3.22/show-coords -rHTl a-b.filter > a-b.coords
rm a-b.delta a-b.filter
perl /ifs4/BC_PUB/biosoft/pipe/bc_mg/BAC_Denovo/BAC_pipeline_1.1.1/Assembly/Assembly_V1.0/lib/synteny/aline_get.pl -series a-b.coords -tlen /hwfssz4/BC_COM_P0/F18FTSECKF1389/ASPjisD/ClustersPloter/example/tmp/filter/Target.len -qlen /hwfssz4/BC_COM_P0/F18FTSECKF1389/ASPjisD/ClustersPloter/example/tmp/filter/Query.len --ident 85 > /hwfssz4/BC_COM_P0/F18FTSECKF1389/ASPjisD/ClustersPloter/example/tmp/filter/Query.len.sort 2> /hwfssz4/BC_COM_P0/F18FTSECKF1389/ASPjisD/ClustersPloter/example/tmp/filter/Target.len.sort
echo finish-a-b > a-b.fin

cd /hwfssz4/BC_COM_P0/F18FTSECKF1389/ASPjisD/ClustersPloter/example/tmp/synteny
/ifs4/BC_PUB/biosoft/pipeline/Package/MUMmer-3.22/nucmer /hwfssz4/BC_COM_P0/F18FTSECKF1389/ASPjisD/ClustersPloter/example/tmp/filter/Target.200.fa /hwfssz4/BC_COM_P0/F18FTSECKF1389/ASPjisD/ClustersPloter/example/tmp/filter/Query2.200.fa -prefix a-c 2>/dev/null
/ifs4/BC_PUB/biosoft/pipeline/Package/MUMmer-3.22/delta-filter -1 a-c.delta > a-c.filter
/ifs4/BC_PUB/biosoft/pipeline/Package/MUMmer-3.22/show-coords -rHTl a-c.filter > a-c.coords
rm a-c.delta a-c.filter
perl -e 'until(-s "a-b.fin"){sleep(5);}'
perl /ifs4/BC_PUB/biosoft/pipe/bc_mg/BAC_Denovo/BAC_pipeline_1.1.1/Assembly/Assembly_V1.0/lib/synteny/aline_get.pl -series -fixing a-c.coords -tlen /hwfssz4/BC_COM_P0/F18FTSECKF1389/ASPjisD/ClustersPloter/example/tmp/filter/Target.len.sort -qlen /hwfssz4/BC_COM_P0/F18FTSECKF1389/ASPjisD/ClustersPloter/example/tmp/filter/Query2.len --ident 85 > /hwfssz4/BC_COM_P0/F18FTSECKF1389/ASPjisD/ClustersPloter/example/tmp/filter/Query2.len.sort

cd /hwfssz4/BC_COM_P0/F18FTSECKF1389/ASPjisD/ClustersPloter/example/tmp/synteny
/ifs4/BC_PUB/biosoft/pipeline/Package/blast-2.2.26/bin/blastall -p blastn -i /hwfssz4/BC_COM_P0/F18FTSECKF1389/ASPjisD/ClustersPloter/example/tmp/filter/Query.200.fa -d /hwfssz4/BC_COM_P0/F18FTSECKF1389/ASPjisD/ClustersPloter/example/tmp/filter/Target.200.fa -e 1e-5 -F F -m 8 | awk '($3>=85 && $4>=100)' > a-b.blastn.m8 
perl /ifs4/BC_PUB/biosoft/pipe/bc_mg/BAC_Denovo/BAC_pipeline_1.1.1/Assembly/Assembly_V1.0/lib/synteny/syncover /hwfssz4/BC_COM_P0/F18FTSECKF1389/ASPjisD/ClustersPloter/example/tmp/filter/Target.len /hwfssz4/BC_COM_P0/F18FTSECKF1389/ASPjisD/ClustersPloter/example/tmp/filter/Query.len a-b.blastn.m8 -rank '1,8,9;0,6,7' > a-b.cover.stat
perl /ifs4/BC_PUB/biosoft/pipe/bc_mg/BAC_Denovo/BAC_pipeline_1.1.1/Assembly/Assembly_V1.0/lib/synteny/stat_coverage.pl -size /hwfssz4/BC_COM_P0/F18FTSECKF1389/ASPjisD/ClustersPloter/example/tmp/filter/Target.len a-b.blastn.m8 -ranks '1,8,9' --len 50 -noann > a-b.cover.lst

cd /hwfssz4/BC_COM_P0/F18FTSECKF1389/ASPjisD/ClustersPloter/example/tmp/synteny
/ifs4/BC_PUB/biosoft/pipeline/Package/blast-2.2.26/bin/blastall -p blastn -i /hwfssz4/BC_COM_P0/F18FTSECKF1389/ASPjisD/ClustersPloter/example/tmp/filter/Query2.200.fa -d /hwfssz4/BC_COM_P0/F18FTSECKF1389/ASPjisD/ClustersPloter/example/tmp/filter/Target.200.fa -e 1e-5 -F F -m 8 | awk '($3>=85 && $4>=100)' > a-c.blastn.m8
perl /ifs4/BC_PUB/biosoft/pipe/bc_mg/BAC_Denovo/BAC_pipeline_1.1.1/Assembly/Assembly_V1.0/lib/synteny/syncover /hwfssz4/BC_COM_P0/F18FTSECKF1389/ASPjisD/ClustersPloter/example/tmp/filter/Target.len /hwfssz4/BC_COM_P0/F18FTSECKF1389/ASPjisD/ClustersPloter/example/tmp/filter/Query2.len a-c.blastn.m8 -rank '1,8,9;0,6,7' > a-c.cover.stat

