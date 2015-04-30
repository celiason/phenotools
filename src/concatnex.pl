#!/usr/bin/perl -w
use diagnostics;
#Brian O'Meara, 2007
#Edited by Chad Eliason, 2014
#concatenatenexusfiles.pl version 1.01
#http://www.brianomeara.info
#Licensed under the GNU Public License Version 2

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

sub concat{
	$labelstates=0;
	print "everything has to be unix-delimited, non-interleaved\n";
	print "file containing list of taxa: ";
	$infile=<>;
	chomp $infile;
	$ntax=0;
	open(IN,"$infile") or die "could not open the file";
	# i think this looks for the species names and collects that list, as well as the data matrix
	while (<IN>) {
		chomp $_;
		push(@taxonarray,"$_");		# collects taxa names
		push(@outputarray,"");		# assign blanks to output array
		$ntax=$ntax+1;
	}
	close IN;
	open(OUT,">concatenated.txt");
	$charsets="begin sets;";
	$charactercount=0;
	# get list of nexus files
	# @nexusfilearray=`ls *.nex`;
	@nexusfilearray = glob('*.nex');
	foreach $nexusfile (@nexusfilearray) {  # the @ denotes an array
		$nexusfile=~m/(.*)\.nex/i;
		$genename=$1;  # gene name given by filename
		$currentcharnumberstring=`grep -i nchar $nexusfile`;  # get number of characters for each file
		$currentcharnumberstring=~m/nchar\s*\=\s*(\d+)/i; # grep character number
		$currentcharnumber=$1;  # use saved number in `()` above
		$rawcodonposset = `grep -i codonposset $nexusfile`;
		$endchar=$charactercount+$currentcharnumber;
		$startchar=$charactercount+1;
		$charsets="$charsets\n"."charset $genename = $startchar - $endchar;";
		open(FH,$nexusfile);
		# chomp(@lines=<FH>);
		# i'm pretty sure the <> around FH indicates to read line by line
		while (<FH>) {
			# collect character labels, looks for this pattern: [1.] ... or [1] ..., store
			# character name and number:
			if (/^.+\[(\d{1,4}|[A-Z]+)[.]*\]\s(.+)/) {
				my ($charnum,$charname) = $_ =~ /^.+\[(\d{1,4})[.]*\]\s'(.+)'/;
				# my ($filename) = $line =~ /(\w+)\.nex/;
				my $newname = "'$charname [$nexusfile trait $charnum]'";
				push(@chararray,$newname);
			}
			# collect state labels:
			if (/STATELABELS/../;$/) {
				if (/^[\t\s]+(\d{1,4})[.]*\s(.+)/){
					my ($statenum,$statename) = $_ =~ /^[\t\s]+(\d{1,4})[.]*\s(.+)/;
					my $newname = "$statename [$nexusfile trait $statenum]";
					push(@statearray,$newname);
				}
			}
		}

# perl -ne 'print if /STATELABELS/ .. /;$/' returns text between markers

# this didn't work:
		# 	if ($line =~ /^.+\d{1,4}\s.+/) {
		# 		my ($statenum,$statename) = $line =~ /^.+(\d{1,4})\s(.+)/;
		# 		my $newname = "$statename [$nexusfile trait $statenum]";
		# 		push(@statearray,$newname);
		# 	}
	# foreach (@lines) {
	# 	if (/TAXLABELS/ ... /\;/) {
	# 		# strip whitespace
	# 		for ($_) {
	# 			s/^\s+//;
	# 			s/\s+$//;
	# 		}
	# 		# add to taxa label array
	# 		push(@taxlabels,$_);
	# 	}
	# 	# get number of taxa
	# 	if (/NTAX[\s]*=[\s]*(\d+)\s*\;/i) {
	# 		$ntax = $1;
	# 	}
	# }
		if ($rawcodonposset=~m/N\:\s*([\d\-\s]+)/) {
			$charsets="$charsets\n"."charset $genename"."intron = ";
			$rawN=$1;
			print "\n$rawcodonposset\t";
			@Narray=split(/\D+/,$rawN);
			foreach $Nelement (@Narray) {
				$newelement=$Nelement+$charactercount;
				$charsets="$charsets"." $newelement";
				print "$newelement\t";
			}
			$charsets="$charsets".";";
		}
		if ($rawcodonposset=~m/1\:\s*([\d\-\s]+)/) {
			$charsets="$charsets\n"."charset $genename"."pos1 = ";
			$raw1=$1;
			@Narray=split(/\D+/,$raw1);
			foreach $Nelement (@Narray) {
				$newelement=$Nelement+$charactercount;
				$charsets="$charsets"." $newelement";
			}
			$charsets="$charsets".";";
		}
		if ($rawcodonposset=~m/2\:\s*([\d\-\s]+)/) {
			$charsets="$charsets\n"."charset $genename"."pos2 = ";
			$raw2=$1;
			@Narray=split(/\D+/,$raw2);
			foreach $Nelement (@Narray) {
				$newelement=$Nelement+$charactercount;
				$charsets="$charsets"." $newelement";
			}
			$charsets="$charsets".";";
		}
		if ($rawcodonposset=~m/3\:\s*([\d\-\s]+)/) {
			$charsets="$charsets\n"."charset $genename"."pos3 = ";
			$raw3=$1;
			@Narray=split(/\D+/,$raw3);
			foreach $Nelement (@Narray) {
				$newelement=$Nelement+$charactercount;
				$charsets="$charsets"." $newelement";
			}
			$charsets="$charsets".";";
		}
		# make data matrix
		# sees how many times taxon listed in nexus file i think
		# this is problematic though in cases where taxon names appear in comments or state descriptions
		# TODO: fix this
		# CE: added \\b so script won't include subspecies (e.g., Homo_sapiens would match Homo_sapiens and Homo_sapiens2)
		# for loop looking when taxon number less than number of taxa within taxonarray (list of all taxa)
		for ($taxonnumber=0;$taxonnumber<scalar(@taxonarray);$taxonnumber++) {
			# see if taxon only listed zero, one times in nexus file
			if (`grep -c '$taxonarray[$taxonnumber]\\b' $nexusfile`<2) {
				for ($i=0;$i<$currentcharnumber;$i++) {
					# puts a question mark if no trait data
					$outputarray[$taxonnumber]="$outputarray[$taxonnumber]"."?";
				}
			}
			elsif (`grep -c '$taxonarray[$taxonnumber]\\b' $nexusfile`==2) {  
				$rawmat=`grep '$taxonarray[$taxonnumber]\\b' $nexusfile`;
				# strip new lines
				$rawmat=~s/\n//g;
				# regular expression that stores everything that isn't a space `(\S+) for use in next line, as `$1`
				if ($rawmat=~m/$taxonarray[$taxonnumber].*$taxonarray[$taxonnumber]\s+(\S+)/) { 
					# `.` is the concatenation operator
					$outputarray[$taxonnumber]="$outputarray[$taxonnumber]"."$1";
				}
				else {
					if($rawmat=~m/$taxonarray[$taxonnumber].*$taxonarray[$taxonnumber](.+)/) {
						print "$taxonarray[$taxonnumber] has "."$1";
					}
					else {
					print "$rawmat\n\n";
					}
				}
			}
			else {
				print "$taxonarray[$taxonnumber] $nexusfile\n";
			}
		}
		$charactercount=$charactercount+$currentcharnumber;
	}

	$charsets="$charsets\nEND;\n";
	print OUT "#NEXUS\nBEGIN TAXA;\nDIMENSIONS NTAX=$ntax;\nTAXLABELS\n";
	$taxonlist=join("\n",@taxonarray);
	print OUT "$taxonlist\n;\nend;\n";
	print OUT "BEGIN CHARACTERS;\nDIMENSIONS NCHAR=$charactercount;\nFORMAT DATATYPE=STANDARD GAP=- MISSING=?;\n";
	# character labeling
	print OUT "CHARLABELS\n";
	for ($charnumber=0;$charnumber<scalar(@chararray);$charnumber++) {
		my $newnum = $charnumber+1;
		print OUT "[$newnum] "."$chararray[$charnumber]\n";
	}
	# chararray is where the character descriptions are
	# my $length = @nexusfilearray;
	# print OUT "there are $length files in the nexus file array";
	# output data matrix
	if ($labelstates==1) {
		print OUT ";\nSTATELABELS\n";
		for ($statenumber=0;$statenumber<scalar(@statearray);$statenumber++) {
			my $newnum = $statenumber+1;
			print OUT "$newnum "."$statearray[$statenumber]\n";
		}
	}
	print OUT ";\nmatrix\n";
	for ($taxonnumber=0;$taxonnumber<scalar(@taxonarray);$taxonnumber++) {
		print OUT "$taxonarray[$taxonnumber]  $outputarray[$taxonnumber]\n";
	}
	print OUT ";\nend;\n$charsets\n";
}







##########################

sub test {
	print "everything has to be unix-delimited, non-interleaved\n";
	print "file containing list of taxa: ";
	$infile=<>;
	chomp $infile;
	$ntax=0;
	open(IN,"$infile") or die "could not open the file";
	# i think this looks for the species names and collects that list, as well as the data matrix
	while (<IN>) {
		chomp $_;
		push(@taxonarray,"$_");		# collects taxa names
		push(@outputarray,"");		# somehow gets numbers after taxon name?
		$ntax=$ntax+1;
	}
	close IN;
	$charsets="begin sets;";
	$charactercount=0;
	$charname="";
	$charnum="";
	@nexusfilearray=`ls *.nex`;
	foreach $nexusfile (@nexusfilearray) {  # the @ denotes an array
		$nexusfile=~m/(.*)\.nex/i;
		$genename=$1;  # gene name given by filename
		$currentcharnumberstring=`grep -i nchar $nexusfile`;  # get number of characters for each file
		$currentcharnumberstring=~m/nchar\s*\=\s*(\d+)/i; # grep character number
		$currentcharnumber=$1;  # use saved number in `()` above
		$rawcodonposset = `grep -i codonposset $nexusfile`;
		$endchar=$charactercount+$currentcharnumber;
		$startchar=$charactercount+1;
		$charsets="$charsets\n"."charset $genename = $startchar - $endchar;";
		open(FH,$nexusfile);
		while (my $line = <FH>) {
			if ($line =~ /^.+\[(\d{1,3}).*\]\s(.+)/) {
				my ($charnum,$charname) = $line =~ /^.+\[(\d{1,3}).*\]\s(.+)/;
				print "$charname, [original trait number was $charnum]\n";	
			}
		}
	}
}

# function call
concat()
