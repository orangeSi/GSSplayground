my $s=1.1;
my $e=1.2;

print "s is $s, e is $e\n";
($s,$e)=&get_real_coordinate($s,$e);
print "s is $s, e is $e\n";

sub get_real_coordinate(){
	my ($s,$e)=@_;
	my $s_precision=length(($s =~ /\.(.*)/)[0]);
	my $e_precision=length(($e =~ /\.(.*)/)[0]);
	print "s_precision $s_precision , e_precision $e_precision\n";
	die "error: precision of $s and $e are not equal\n" if($s_precision != $e_precision);
	my $unit;
	if($s == $e){
		$unit=1/(10**$s_precision);
		$s=$s;
		$e=$s+$unit;
		print "1unit is $unit\n";
	}elsif($s<=$e){
		$unit=1/(10**$s_precision);
		$s=$s;
		$e=$e+$unit;
		print "2unit is $unit\n";
	}else{
		die "error:$s should <= $e, but not in fact\n"
	}
	return ($s,$e);
}
