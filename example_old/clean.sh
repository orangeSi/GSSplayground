set -vx
ls *notitle*svg|while read svg
do
	clean=`echo $svg|sed 's/notitle/cleaned/'`
	svgcleaner $svg $clean 2>/dev/null && echo get $cleaned.svg && rm $svg

done
