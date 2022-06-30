#!/usr/bin/perl -w

use strict;
use Statistics::Descriptive;

my $input = "";
my $output = "yn00_partitions.log";
#my $split = 0;
my $precision = 3;

if ($#ARGV == -1) {
	print "\n";
	print "This script parses output from codeml runmode=-2, extracts the dN and dS values and calculates the mean and SD.\n";
	print "\nSyntax:\n";
	print " -f X		use file X as input\n";
	print " -out X		output to file X\n";
#	print " -split X	use sequence number X as first in the new group\n";	
	print " -precision X	set precision to X decimal places (3 by default)\n";
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
#	if ($ARGV[$argnum] eq "-split" and $ARGV[$argnum+1]) {
#		$split = $ARGV[$argnum+1];
#	}
	if ($ARGV[$argnum] eq "-precision" and $ARGV[$argnum+1]) {
		$precision = $ARGV[$argnum+1];
	}
}

if ($input) {
	open (CODEML_FILE, $input);
	open (LOG, ">$output");
	
	print "\nLoading file $input...\n";
	print LOG "\nLoading file $input...\n";
	
	my $dnds_block = 0;

	my $bg_dS = Statistics::Descriptive::Full->new();
	my $bg_dN = Statistics::Descriptive::Full->new();

#	my $fg_dS = Statistics::Descriptive::Full->new();
#	my $fg_dN = Statistics::Descriptive::Full->new();
	while (<CODEML_FILE>) {
		chomp;
		if ($_ eq "") {
			next;
		}
		if ($_ =~ "pairwise comparison, codon frequencies") {
			$dnds_block = 1;
			next;
		}
#		if ($_ eq "(C) LWL85, LPB93 & LWLm methods") {
#			$yn_block = 0;
#			last;
#		}
		if ($dnds_block == 1) {
			my @fields = split(' ');
			if ($fields[0] eq "dN" and $fields[1] eq "=" and $fields[12] eq "1)") {
#				if ($split != 0 and $fields[0] >= $split and $fields[1] >= $split) {
#					$fg_dS->add_data($fields[10]); 
#					$fg_dN->add_data($fields[7]); 
#				}
#				elsif ($split == 0 or ($fields[0] < $split and $fields[1] < $split)) {
					$bg_dS->add_data($fields[7]);
					$bg_dN->add_data($fields[2]); 				
#				}
			}
		}
	}
	my $fmt = "%.".$precision."f";
	
	my $bg_dS_mean = sprintf($fmt, $bg_dS->mean());
	my $bg_dS_sd = sprintf($fmt, $bg_dS->standard_deviation()); 

	my $bg_dN_mean = sprintf($fmt, $bg_dN->mean());
	my $bg_dN_sd = sprintf($fmt, $bg_dN->standard_deviation()); 

	my $bg_count = $bg_dS->count();
	print "\nBackground group mean dS: $bg_dS_mean +/- $bg_dS_sd ($bg_count items)\n";
	print "Background group mean dN: $bg_dN_mean +/- $bg_dN_sd ($bg_count items)\n";	
	
	print LOG "\nBackground group mean dS: $bg_dS_mean +/- $bg_dS_sd ($bg_count items)\n";
	print LOG "Background group mean dN: $bg_dN_mean +/- $bg_dN_sd ($bg_count items)\n";	
	
#	if ($split > 0) {
#		my $fg_dS_mean = sprintf($fmt, $fg_dS->mean());
#		my $fg_dS_sd = sprintf($fmt, $fg_dS->standard_deviation());
#		
#		my $fg_dN_mean = sprintf($fmt, $fg_dN->mean());
#		my $fg_dN_sd = sprintf($fmt, $fg_dN->standard_deviation()); 
#		 
#		my $fg_count = $fg_dS->count();
#		print "\nForeground group mean dS: $fg_dS_mean +/- $fg_dS_sd ($fg_count items)\n";
#		print "Foreground group mean dN: $fg_dN_mean +/- $fg_dN_sd ($fg_count items)\n";
#		
#		print LOG "\nForeground group mean dS: $fg_dS_mean +/- $fg_dS_sd ($fg_count items)\n";
#		print LOG "Foreground group mean dN: $fg_dN_mean +/- $fg_dN_sd ($fg_count items)\n";
#	}
	print "\nFinished reading, exiting.\n";
	print LOG "\nFinished reading, exiting.\n";			
	
	close LOG;
	close CODEML_FILE;

}
else {
	print "\nNo file specified. Use the \"-f\" param to pass a filename.\n";
}

exit;

