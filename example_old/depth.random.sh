sample="s3"
scf="s3"
avg_dep=50
len=10000
out=$sample.$scf.depth.txt
echo -e "#sample\tscf_id\tpos\tdepth" > $out
for i in $(seq 1  $len)
do
    depth=$(( ( RANDOM % $avg_dep )  + 1 ))
    echo -e "$sample\t$scf\t$i\t$depth" >>$out
done
echo out is $out
