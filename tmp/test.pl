use myth qw(sub1 sub2 sub3 sub4 sub5);

my %a;
$a{"1"}=2;
print "first $a{1}=2\n\n";

my @arr1=(1,2,3,4);
my @arr2=(5,6,7,8);

&sub1("dd", \%a);


my ($a, $a2) = &sub2(%a);
@a2=@$a2;
#die "error: a2 is @a2\n";

&sub4(\@arr1,\@arr2);
print "dd is $a{2}\n";
print "done\n";

my $tmp=1;
&sub3($tmp, \%a);
print "tmp is $tmp\n";

@arr1=(1,2,3,4);
@arr2=(11,21,31,14);
&sub5("tt", @arr1, @arr2);
