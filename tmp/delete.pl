my %hash = (foo => 11, bar => 22, baz => 33);
print keys(%hash);
print "\n";
delete $hash{foo}; 
print keys(%hash);
