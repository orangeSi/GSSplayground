
for blast in NC_017659.1.fa.to.NC_017659.1.fa.blast.m6 NC_017659.1.indel.fa.to.NC_017659.1.fa.blast.m6
do
	echo $blast
	echo -e "#query\ttarget\tquery_start\tquery_end\ttarget_start\ttarget_end\tstrand\tindentity(just_for_display)\talignment_length(just_for_display)\tpopup" >$blast.common
	cat $blast|awk -F '\t' '{t="+";if($9>$10){t="-"};print $1,$2,$7,$8,$9,$10,t,$3,$4,""}'|tr ' ' '\t' >>$blast.common
	echo $blast.common
	echo
done
