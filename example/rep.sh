set -vex
for i in $(seq 1 1 8)
do
#	sh ../ClustersPloter.sh tracks.$i.list out$i . main.$i.conf && cp out$i.html out$i.rep1.html && cp out$i.svg  out$i.rep1.svg
#	sh ../ClustersPloter.sh tracks.$i.list out$i . main.$i.conf && cp out$i.html out$i.rep2.html && cp out$i.svg  out$i.rep2.svg
	diff out$i.rep1.html out$i.rep2.html|wc -l
done

