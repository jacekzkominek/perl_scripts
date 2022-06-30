#!/usr/bin/perl -w
#SAMPLE ALICAT COMMAND:
#alicat.pl -s -o SSB1_scer.fas -r SSB1:/media/sda4/work/zestaw/polymorphism_data/Scerevisiae_Sparadoxus/misc/cere/ref/genome.gff /media/sda4/work/zestaw/polymorphism_data/Scerevisiae_Sparadoxus/cere/match/*/imputed.gz

use strict;

my $input = "";
my $output = "outfile.fas";

if ($#ARGV == -1) {
	print "\n";
	print "This script parses output from the alicat.pl script from SGRP polymorphism data.\n";

	print "\nSyntax:\n";
	print " -f X		use file X as input\n";
	print " -out X		output to file X\n";
	print "\n";
	exit;
}

foreach my $argnum (0 .. $#ARGV) {
	if ($ARGV[$argnum] eq "-f" and $ARGV[$argnum+1]) {
		$input = $ARGV[$argnum+1];
	}
	if ($ARGV[$argnum] eq "-out" and $ARGV[$argnum+1]) {
		$output = $ARGV[$argnum+1];
	}
}

if ($input) {
	print "\nLoading file $input...";
	open (ALICAT_FILE, $input) or die "Cannot open $input";
	print "\nParsing strain headers.";

	my @heads;
	my @seqs;
	while (<ALICAT_FILE>) {
		chomp;
		if ((substr $_,0,1) eq ">") {
			@heads = split(' ');
			shift @heads;
		}
	}
	seek ALICAT_FILE, 0, 0;
	for (my $i1=0;$i1 < @heads;$i1++) {
		$seqs[$i1] = "";
	}
	print "\nParsing strain sequences.";
	while (<ALICAT_FILE>) {
		my @fields = split (//);
		if (@fields >= 17 and $fields[16] eq '|') {
			for (my $i1=0;$i1 < @heads;$i1++) {
				if (@fields >= 17+$i1) {
					$seqs[$i1] = $seqs[$i1].$fields[17+$i1];
				}
			}
		}
	}	
	close ALICAT_FILE;
	print "\nWriting to file $output.\n\n";
	open (OUT, ">$output");
	for (my $i1=0;$i1 < @heads;$i1++) {
		print OUT ">".$heads[$i1]."\n";
		print OUT uc($seqs[$i1])."\n";
	}	
}
else {
	print "\nNo file specified. Use the \"-f\" param to pass a filename.\n";
}

exit;

