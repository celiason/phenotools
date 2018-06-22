#!/usr/bin/perl -w

# use strict;
# use warnings;

my $filename = "toy1.nex";
#open(FH, $filename);
#my @lines = <FH>;
#close(FH);

#my $count = scalar @lines;  # scalar counts number of lines in the file
#print("There are $count lines in $filename\n");

sub printie
{
	open(FH, $filename);
	while(my $line = <FH>) {
		if($line=~/^.\[\d{1,2}\](.*)/) {
			print("name is $1\n")
			}
	}
}


#my $charlabels =~m/\[\d{1,2}\]+/;

#print(m/\[\d{1,2}\]+/);

#print("$charlabels");

main(@ARGV);

# printie(@ARGV);


sub main
{
	my $num = 0;
	# my @array = (1, "two", 3, 4);
	open(FH, $filename);
	while (<FH>) {
		push(@array,/^.\[\d{1,2}\](.*)/);		# collects character labels
	}
	# while(my $line = <FH>) {
	# 	if($line=~/^.\[\d{1,2}\](.*)/) {
	# 		print("name is $1\n")
	# message("there are" . scalar @array . "elements in the array");
	# print(join('\n', @array));
	# print("this is my array - @array\n");
	
	$charlist=join("\n",@array);
	print("$charlist\n");
	# for ($num=0;$num<scalar(@array);$num++)
	# { 
	# 	print "$array[$num]\n";
	# }
}



