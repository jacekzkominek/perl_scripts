#!/usr/bin/perl -w

use strict;
use Bio::AlignIO;
use File::Find::Object::Rule;

my $in_type = "Nexus";
my $in_ext = "nex";
my $out_type = "Fasta";
my $out_ext = "fas";
    
if ($#ARGV == -1) {
	print "\n";
	print "This script converts a sequence file from $in_type (*.$in_ext) to $out_type (*.$out_ext) format.\n";
	print "\nSyntax:\n";
	print " -f X		use file X as input or 'all' to convert all files in the current dir\n";
	print " -r 		recursively parse all dirs below the current (only usable with '-f all')\n";
	print "\n";
	exit;
}
my $input = "";
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
		my $ffor = File::Find::Object::Rule->file()->name("*.$in_ext");
		my $str = " (and below)";
		if ($recurse == 0) {
			$ffor->maxdepth(1);
			$str = "";
		}
		print "\nConverting all $in_type files from current directory$str to $out_type format...\n";	

		my @filelist = $ffor->in(".");
		foreach my $file (@filelist) {
			convert_file($file);
		}
	}
	elsif (-e $input) {
		print "\nConverting $input to $out_type format...\n";	
		convert_file($input);
	}
}
else {
	print "\nNo file specified or file does not exists. Pass a filename after the \"-f\" parameter (or 'all' to use all files in the current dir)\n";
}
exit;


sub convert_file {
	my $file = $_[0];
	my $seq_in = Bio::AlignIO->new(-file => $file, -format => "$in_type");
	my $file_short = substr($file,0,length($file)-4);
	my $seq_out = Bio::AlignIO->new(-file => ">$file_short.$out_ext", -format => "$out_type", -show_symbols => 0, -show_endblock => 0);	
	while (my $align = $seq_in->next_aln() ) {
		$seq_out->write_aln($align);
    	}
}

