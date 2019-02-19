set -vex
pdi=150
ls out*notitle.svg|while read line
do
	png=`echo $line|sed 's/svg$/png/'`
	convert -density $pdi $line $png
done
