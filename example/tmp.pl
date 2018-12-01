

my %fs;
$fs{"1"}=2;

&A();

sub A(){
    print "A\n";
    for my $k(keys %fs){
        print "$k->$fs{$k}\n";
    }
    &B();
}
sub B(){
    print "B\n";
    for my $k(keys %fs){
        print "$k->$fs{$k}\n";
    }
}
