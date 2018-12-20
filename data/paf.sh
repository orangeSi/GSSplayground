minimap2 s2.seq s1.seq -cx asm20 |gzip >s1.mapto.s2.paf.gz
minimap2 s3.seq s2.seq -cx asm20 |gzip >s2.mapto.s3.paf.gz
minimap2 s4.seq s2.seq -cx asm20 |gzip >s2.mapto.s4.paf.gz
