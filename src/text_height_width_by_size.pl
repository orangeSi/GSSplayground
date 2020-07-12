use lib "./";
use JSON;
use List::Util qw(max min);
use Imager::Font;

my $out="text_height_width_by_size.json";
my %height_width;
my $text = "W";
my $ttf="Times_New_Roman.ttf";
die "error: ttf $ttf not exists" if(!-f $ttf);
print "using $ttf\n";
my $font = Imager::Font->new(file => $ttf);

for my $int(0..100){
	for my $step(0..10){
		$fontsize=$int+$step/10;
		print "text is $text, fontsize is $fontsize, ttf is $ttf\n";
		my $bbox = $font->bounding_box(string=>"$text", size=>$fontsize);
		$height_width{$fontsize}{height}{min} = min($bbox->font_height, $bbox->text_height);
		$height_width{$fontsize}{height}{max} = max($bbox->font_height, $bbox->text_height);
		$height_width{$fontsize}{height}{"avg"} = ($bbox->font_height+$bbox->text_height)/2;

		$height_width{$fontsize}{width} = $bbox->total_width;

#if($type eq "height"){
#	return min($bbox->font_height, $bbox->text_height) if($value eq "min");
#	return max($bbox->font_height, $bbox->text_height) if($value eq "max");
#	return ($bbox->font_height+$bbox->text_height)/2;
#}elsif($type eq "width"){
#	return $bbox->total_width;
#}else{
#	die "error: not support $type for check_font_size_by_estimate， only height or width"
#}
	}
}
open OUT,">$out";
print OUT encode_json(\%height_width);
close OUT;
print("write to $out\n");


open IN,"$out";
my $myHashEncoded=<IN>;chomp $myHashEncoded;
my $myHashRefDecoded = decode_json($myHashEncoded);
my %myHashEncoded = %$myHashRefDecoded;

my $fs=99;
for my $key(keys %{$myHashEncoded{$fs}{height}}){
	print "$key -> $myHashEncoded{$fs}{height}{\"$key\"}->$myHashEncoded{$fs}{width}\n";

}

my $fs="99";
for my $key(keys %{$myHashEncoded{$fs}{height}}){
	print "$key -> $myHashEncoded{$fs}{height}{\"$key\"}->$myHashEncoded{$fs}{width}\n";

}
my $fs="99.9";
for my $key(keys %{$myHashEncoded{$fs}{height}}){
	print "$key -> $myHashEncoded{$fs}{height}{\"$key\"}->$myHashEncoded{$fs}{width}\n";

}
