use Storable;
#my %table;
#$table{'1'}=2;
#store \%table, 'file';
my $hashref = retrieve('file');
print keys(%$hashref);
