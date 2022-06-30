#!/usr/bin/perl -w

use strict;
use IO::String;
use Bio::SearchIO;
use Cwd;
use File::Find::Object::Rule;

my $input="";
my $target_dir = getcwd;

if ($#ARGV == -1) {
	print "\n";
	print "This script summarizes output from blast_parse.pl.\n";
	print "\nSyntax:\n";
	print " -f X		use file X as input or 'all' to parse all files in the current dir\n";
	print " -dir X		save output in dir X\n";	
	print "\n";
	exit;
}

foreach my $argnum (0 .. $#ARGV) {
	if ($ARGV[$argnum] =~ "-f") {
		$input = $ARGV[$argnum+1];
	}
	if ($ARGV[$argnum] =~ "-dir") {
		$target_dir = $ARGV[$argnum+1];
	}
}
if ($input) {
	print "\nParsing BLAST logs...\n";
	open(LOG, ">$target_dir\/BLAST_log.csv"); 
	print LOG "Hit_name,Score,E-value\n";
	if ($input eq "all") {
		my $ffor = File::Find::Object::Rule->file()->name(/blast.+\.txt/);
		$ffor->maxdepth(1);
		my @filelist = $ffor->in(".");
		foreach my $file (@filelist) {
			parse_logfile($file);
		}
	}
	elsif (-e $input) {
		parse_logfile($input);
	}
	close(LOG);
}
else {
	print "\nNo file specified or file does not exists. Pass a filename after the \"-f\" parameter (or 'all' to use all files in the current dir)\n";
}
exit;

sub parse_logfile {
	my $file = $_[0];
	print "Parsing - $file...";
	my $report = Bio::SearchIO->new(-file => $file);
	my $int = 0;
	my $result = $report->next_result();
	my $db_name = $result->database_name();
	print LOG "$db_name\n";
	while ($int < $result->num_hits and $int < 5) {
		my $hit = $result->next_hit();
#		$hit_name = $hit->name().$hit->description();
		my $hit_name = $hit->name();
		my $hit_score = $hit->score();
		my $hit_e = $hit->significance();
		print LOG "$hit_name,$hit_score,$hit_e\n";
		$int++;
	}
	print "finished\n";
}

