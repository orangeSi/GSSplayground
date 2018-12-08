
my ($cg)=@_;
die "error: cg is $cg\n" if($cg!~ /^\d+[^\d]+$/);
$cg=~ /(\d+)([^\d]+)/;
my $cg_len=$1;
my $cg_type=$2;
my $color_height_cs="M:green:0.5:1bp,I:red:0.7:6bp,D:black:0.7:5bp,N:blue:0.2:1bp,S:blue:0.2:1bp,H:blue:0.2:1bp,P:blue:0.2:1bp,X:grey:0.6:1bp,reverse:#D2691E:0.8:6bp";
my $color_height_cs_usage="M:green:0.5:1bp,I:red:0.7:6bp,D:black:0.7:5bp,N:blue:0.2:1bp,S:blue:0.2:1bp,H:blue:0.2:1bp,P:blue:0.2:1bp,X:grey:0.6:1bp,reverse:#D2691E:0.8:6bp";
my (%colors_height);
my ($cs_color, $cs_height);
$color_height_cs=~ s/\s//g;
my @color_height_cses=split(/,/, $color_height_cs);
for my $ch(@color_height_cses){
	my @arr=split(/:/, $ch);
	die "error:$ch of $color_height_cs format is wrong, should like $color_height_cs_usage\n" if(@arr!=4 || $arr[-1]!~ /bp$/);
	my ($cg,$color,$height,$limit_len)=@arr;
	$limit_len=~ s/bp//g;
	$colors_height{$cg}{color}=$color;
	$colors_height{$cg}{height}=$height;
	$colors_height{$cg}{limit_len}=$limit_len;
}

my $feature_color;
my $feature_height;
my $cg_before="";
my @cigars=$cg=~ /([\d\+]+[^\d^\+])/g;
$cigars_len=scalar(@cigars);
for my $i(0..$cigars_len){
	$i=~ /(\d+)([^\d])/;
	my ($cg_len, $cg_type)= ($1,$2);
	die "error:not support cigar $cg_type for $i\n" if(not exists $colors_height{$cg_type});
	if($cg_len < $colors_height{$cg_type}{limit_len}){
		$cigars[$i]="$cg_len"."M";
	}
}
$cg=join "", @cigars;
while($cg_before ne $cg){
	$cg_before=$cg;
	$cg=~ s/(\d+)M(\d+)M/$1+$2M/g;
}
@cigars=$cg=~ /([\d\+]+[^\d^\+])/g;
$cigars_len=scalar(@cigars);
for my $i(0..$cigars_len){
	$cigars[$i]=~ /([\d\+]+)([^\d^\+])/;
	my ($cg_len, $cg_type)= ($1,$2);
	if($cg_len=~ /\+/){
		$cg_len=eval($cg_len);
		$cigars[$i]="$cg_len$cg_type";
	}
}
$cg=join "", @cigars;
