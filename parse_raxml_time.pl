#!/usr/bin/perl -w

use strict;
use Statistics::Descriptive;
use File::Find::Object::Rule;

my $input = "";

if ($#ARGV == -1) {
	print "\n";
	print "This script parses RAxML log files and outputs bootstrap and total time of the analysis..\n";
	print "\nSyntax:\n";
	print " -f X		use file X as input or 'all' to parse all files in the current dir\n";
	print "\n";
	exit;
}

my $recurse = 0;

foreach my $argnum (0 .. $#ARGV) {
	if ($ARGV[$argnum] eq "-f" and $ARGV[$argnum+1]) {
		$input = $ARGV[$argnum+1];
		if (-e $input) {
			print "\nLoading $input...\n";
		}
		elsif ($input ne "all") {
			$input="";
		}
	}
	if ($ARGV[$argnum] eq "-r") {
		$recurse = 1;
	}
}

if ($input) {
	if ($input eq "all") {
		my $ffor = File::Find::Object::Rule->file()->name("RAxML_info.*");
		$ffor->maxdepth(1);
		my @filelist = $ffor->in(".");
		my $times = Statistics::Descriptive::Full->new();
		foreach my $file (@filelist) {
			parse_raxml_log($file,$times);
		}
		my $times_mean = $times->mean();
		my $times_SD = $times->standard_deviation();
		print "\nMean ratio: $times_mean +/- $times_SD\n";
	}
	elsif (-e $input) {
		parse_raxml_log($input);
	}
}
else {
	print "\nNo file specified or file does not exists. Pass a filename after the \"-f\" parameter (or 'all' to use all files in the current dir)\n";
}

exit;

sub parse_raxml_log {
	my $log_file = $_[0];
	my $time_stat = $_[1];
	print "\nParsing $log_file...";
	open (RAXML_LOG, $log_file);
	my $bs_avg_time = 0;
	my $total_time = -1;
	while (<RAXML_LOG>) {
		chomp;
		if ($_) {
			my @fields = split(' ');
			if ($fields[0] eq "Average") {
				$bs_avg_time = $fields[5];
			}
			if ($fields[0] eq "Overall" and $fields[1] eq "execution") {
				$total_time = $fields[7];		
			}
		}
	}
	my $bs_total_time = $bs_avg_time*100;
	my $ratio = ($bs_total_time)/$total_time;
	$time_stat->add_data($ratio);
	print "\nBootstrap-to-total time ratio: $ratio\n";

	close (RAXML_LOG);
}

