#!/usr/bin/env perl -w
##use strict;
use warnings;

my $t=1;

print "t is $t\n";
&xx();
print "t is $t\n";
sub xx(){
	$t=2;
}
