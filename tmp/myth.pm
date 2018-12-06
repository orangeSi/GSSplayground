package myth;
use strict;
use warnings;

use Exporter qw(import);

our @EXPORT_OK = qw(sub1 sub2 sub3 sub4 sub5);

sub sub1(){
    my ($tmp, $a)=@_;
    &sub3($tmp, $a);
    print "sub1 $_[0]\n";
    print "sub1 again $_[0]\n";
    print "$a->{1}=2\n\n";
    
}
sub sub2(){
    my (%a)=@_;
    print "sub2\n";
    print "$a{1}=2\n\n";
    $a{2}=3;
    my @a2=(1,2,3);
    return (\%a,\@a2);
}

sub sub3(){
    my ($tmp, $a)=@_;
    print "sub3\n";
    print "$a->{1}=2\n\n";
    $tmp=333;
}

sub sub4(){
    my ($arr1, $arr2)=@_;
    print "sub5\n";
    print "arr1 is @{$arr1}; arr2 is @{$arr2}\n";
#    my @arr1=@{$_[0]};
#    my @arr2=@{$_[1]};
#    print "sub5\n";
#    print "arr1 is @arr1; arr2 is @arr2\n";
}

sub sub5(){
    my ($tmp, @a1, @a2)=@_;
    print "sub55\n@a1\n";
    print "sub55\n@a2\n";
}

1;
