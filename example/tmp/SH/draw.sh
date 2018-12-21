### for nucleotide synteny
cd /hwfssz4/BC_COM_P0/F18FTSECKF1389/ASPjisD/ClustersPloter/example/tmp/figure
## draw XOY synteny:
perl /ifs4/BC_PUB/biosoft/pipe/bc_mg/BAC_Denovo/BAC_pipeline_1.1.1/Assembly/Assembly_V1.0/lib/synteny/alig_line.pl -idencut 85 -trade /hwfssz4/BC_COM_P0/F18FTSECKF1389/ASPjisD/ClustersPloter/example/tmp/synteny/a-b.blastn.m8 -rank '1,8,9;0,6,7,12,2' -tlen /hwfssz4/BC_COM_P0/F18FTSECKF1389/ASPjisD/ClustersPloter/example/tmp/filter/Target.len.sort -qlen /hwfssz4/BC_COM_P0/F18FTSECKF1389/ASPjisD/ClustersPloter/example/tmp/filter/Query.len.sort -name -xtitle Target -ytitle Query -xgrid -ygrid -symbol > Target-Query.nuc.xoy.svg
`convert Target-Query.nuc.xoy.svg Target-Query.nuc.xoy.png`
perl /ifs4/BC_PUB/biosoft/pipe/bc_mg/BAC_Denovo/BAC_pipeline_1.1.1/Assembly/Assembly_V1.0/lib/synteny/alig_line.pl -idencut 85 -trade /hwfssz4/BC_COM_P0/F18FTSECKF1389/ASPjisD/ClustersPloter/example/tmp/synteny/a-c.blastn.m8 -rank '1,8,9;0,6,7,12,2' -tlen /hwfssz4/BC_COM_P0/F18FTSECKF1389/ASPjisD/ClustersPloter/example/tmp/filter/Target.len.sort -qlen /hwfssz4/BC_COM_P0/F18FTSECKF1389/ASPjisD/ClustersPloter/example/tmp/filter/Query2.len.sort -name -xtitle Target -ytitle Query2 -xgrid -ygrid -symbol > Target-Query2.nuc.xoy.svg
`convert Target-Query2.nuc.xoy.svg Target-Query2.nuc.xoy.png`
## draw parallel synteny figure:
perl /ifs4/BC_PUB/biosoft/pipe/bc_mg/BAC_Denovo/BAC_pipeline_1.1.1/Assembly/Assembly_V1.0/lib/synteny/aligment_plot.pl -idencut 85 -trade /hwfssz4/BC_COM_P0/F18FTSECKF1389/ASPjisD/ClustersPloter/example/tmp/synteny/a-b.blastn.m8 /hwfssz4/BC_COM_P0/F18FTSECKF1389/ASPjisD/ClustersPloter/example/tmp/synteny/a-c.blastn.m8 -ab_rank '1,8,9;0,6,7,12,2' -ac_rank '1,8,9;0,6,7,12,2' -alen /hwfssz4/BC_COM_P0/F18FTSECKF1389/ASPjisD/ClustersPloter/example/tmp/filter/Target.len.sort -blen /hwfssz4/BC_COM_P0/F18FTSECKF1389/ASPjisD/ClustersPloter/example/tmp/filter/Query.len.sort -clen /hwfssz4/BC_COM_P0/F18FTSECKF1389/ASPjisD/ClustersPloter/example/tmp/filter/Query2.len.sort -grid -name -symbol 'Forward chain[,Reverse chain[,Forw alignment,Reve alignment' -species Target,Query,Query2 > Target-Query-Query2.nuc.parallel.svg
`convert Target-Query-Query2.nuc.parallel.svg Target-Query-Query2.nuc.parallel.png`
