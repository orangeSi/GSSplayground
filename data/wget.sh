id="NC_017659.1"
wget "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nucleotide&id=$id&rettype=fasta" -O $id.fa
wget ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/183/345/GCF_000183345.1_ASM18334v1/GCF_000183345.1_ASM18334v1_genomic.gff.gz -O tmp.gz && gzip -dc tmp.gz|grep $id >$id.gff && rm tmp.gz
