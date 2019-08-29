



echo -e "#track_name\tgff\tgenome_or_length\tchr_id\tstart\tend\tchr_id\tstart\tend\t..." > mytrack.list
echo -e "mychr1\tchr1.gff\tchr1.fa.length\tchr1\t1000\t3000\tchr1\t5000\t8000" >>mytrack.list
echo -e "chr1\t10000" >chr1.fa.length
echo -e "chr1\tfake\tgene\t1300\t1600\t.\t+\t.\tID=gene1;xxx=yyy;" >chr1.gff
echo -e "chr1\tfake\tgene\t1800\t2600\t.\t+\t.\tID=gene2;xxx=yyy;" >>chr1.gff
echo -e "chr1\tfake\tgene\t5700\t6000\t.\t+\t.\tID=gene3;xxx=yyy;" >>chr1.gff


echo -e "feature_keywords = gene, " >mymain.conf
echo -e "feature_setting = myfeature_setting.txt" >>mymain.conf

echo -e "gene1\tfeature_color\tblue" >myfeature_setting.txt
echo -e "gene3\tfeature_color\tblack" >>myfeature_setting.txt
echo -e "gene2\tfeature_color\tgreen" >>myfeature_setting.txt

