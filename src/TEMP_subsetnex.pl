#!/usr/bin/perl -w
use diagnostics;
#Brian O'Meara, 2007
#Edited by Chad Eliason, 2014
#concatenatenexusfiles.pl version 1.01
#http://www.brianomeara.info
#Licensed under the GNU Public License Version 2

#  Script to sort nexus data by character labels (e.g., to have all 'rostrum' traits together)

#  This script should be placed in a folder containing individual datasets for each gene,
#  each ending in ".nex", and a file, not ending in ".nex", which contains a list of
#  all the taxa (each separated by a line feed; no spaces in taxon names). The names
#  of the taxa must match between this file and the gene files, though not every
#  taxon need be sequenced for every gene. Upon executing this script, it will
#  ask you for the name of the file of taxa, then use this and all the *nex files
#  to create a new nexus file containing all the sequences, with character sets
#  automatically created based on the input genes. If the gene files come
#  from Mesquite or MacClade, with codon positions stored, this information
#  is used to create  character sets for each codon position of each gene.
#  Feel free to modify the script as needed.
#  Note that the input files should have unix line breaks and the gene files
#  must not be interleaved.


# TODO add ability to sort

sub test{
	my $ntax;
	my @taxalabels;
	print "nexus file containing data: ";
	chomp ($infile=<>);
	open IN, $infile;
	# my $lines = join '', <IN>;
	chomp (@lines = <IN>);
	# print $lines;
	# my ($match) = $lines =~ /NTAX[\s]*=[\s]*(\d+)\s*\;/ig;
	# my (@taxlabels) = $lines =~ /TAXLABELS(.+)\;/g;
	# my (@taxlabels) = $lines =~ /TAXLABELS[\s\n\w]+/;
	foreach (@lines) {
		if (/TAXLABELS/ ... /\;/) {
			# strip whitespace
			for ($_) {
				s/^\s+//;
				s/\s+$//;
			}
			# add to taxa label array
			push(@taxlabels,$_);
		}
		# get number of taxa
		if (/NTAX[\s]*=[\s]*(\d+)\s*\;/i) {
			$ntax = $1;
		}
	}
	# test if it works:
	# my $x = chomp(@taxlabels);
	# print $x;
	# print "@taxlabels\n";
	print "@taxlabels[0] is the first one\n";
	print "There are $ntax taxa\n";
	close IN;
}

# bertelli_2002.nex

# test();

sub matchnex{
	print "nexus file containing data: ";
	chomp ($infile=<>);
	open(IN,$infile);

	# get data

	# get char labels

	# sort char labels

	# rearrange data based on sorted values of char labels

	# output charlabels, data sorted


	# look for text within character labels
		# if($rawmat=~m/$taxonarray[$taxonnumber].*$taxonarray[$taxonnumber]\s+(\S+)/) { 
					# `.` is the concatenation operator
					# $outputarray[$taxonnumber]="$outputarray[$taxonnumber]"."$1";
				# }

	while (<IN>) {
		# my $x = /.*\w+[\s\t]*(.+)/;
		my $pat = '(\d{1}|\(\d{1,3}\)|\-|\?)';		
		my $x= '1(10)101(01)000001111';
		my @c = $x =~ /$pat/g;
		# my $x = "Hi chad";
		my $count = @c;
		print "-------------\n";
		print "The pattern $pat matches $count times in $x\n";
		print "Here are the matches: @c\n";
		print "Here is match #2: $c[1]\n";
		my @c2 = sort(@c);
		print "Here they are sorted: @c2\n";
		print "-------------\n";
	}
	close(IN);
}

matchnex();



# TODO add ability put close together based on 'fuzzy' match between two names

sub test2{
	# declare variables, arrays, etc.
	my @subset;
	my @traitnum;
	my $text="premax";
	# my $match;
	# load file
	print "nexus file containing data: ";
	$infile=<>;
	chomp $infile;
	open(IN,$infile);
	# look for text within character labels
	while (<IN>) {
		if (/.+\[\d{1,4}[.]*\]\s.+/) {
			if (/$text/) {
				my ($num,$match) = /(\d{1,4})(.+)/;
				push(@subset,"$match");
				push(@traitnum,$num);
			}
		}
	}
	print "@subset\n";
	print "@traitnum\n";
	close(IN);
}
