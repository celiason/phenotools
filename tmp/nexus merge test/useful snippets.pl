#!/usr/bin/perl
# parse_genes.pl
use strict; use warnings;

while (my $line = <>) {
	if ($line =~ /^\s{5}gene/) {
		my ($beg, $end) = $line =~ /(\d+)\.\.(\d+)/;
		$line = <>;
		my ($name) = $line =~ /="(.+)"/;
	 	print "$name $beg $end\n";
	 	}
}