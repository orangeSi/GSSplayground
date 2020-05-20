set -vx
for i in $(seq 1 1 8)
do
	#sh ../ClustersPloter.sh tracks.$i.list out$i . main.$i.conf
	svgcleaner out$i.svg out$i.cleaned.svg 2>/dev/null && echo get out$i.cleaned.svg && rm out$i.notitle.svg
	#cairosvg out$i.cleaned.svg   -o  out$i.cleaned.cairosvg.pdf

done
ls -lhtr *error 
ls -lhtr *html
