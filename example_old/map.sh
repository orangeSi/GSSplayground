/ifshk7/BC_PS/sikaiwei/assembly/miniasm/minimap2/minimap2-2.11_x64-linux/minimap2 -x asm5  s1.seq s2.seq -c >s1.map.to.s2.paf
less s1.map.to.s2.paf |awk -v prefix="fake" '{print $1"\tGenBank\tgene\t"$3"\t"$4"\t.\t.\t.\tID="$1"."prefix"."NR";\n"$6"\tGenBank\tgene\t"$8"\t"$9"\t.\t.\t.\tID="$6"."prefix"."NR";"}' >fake.gff
cat s1.gff fake.gff >s1.add.fake.gff
cat s2.gff fake.gff >s2.add.fake.gff
