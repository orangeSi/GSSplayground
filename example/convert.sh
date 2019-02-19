set -vex
pdi=100
ls out*notitle.svg|while read line
do
	png=`echo $line|sed 's/svg$/png/'`
	convert -density $pdi $line $png
done
