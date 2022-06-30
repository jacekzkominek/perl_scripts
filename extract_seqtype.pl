#!/usr/bin/perl -w

use strict;
use Bio::Seq;
use Bio::SeqIO;
use File::Find::Object::Rule;

if ($#ARGV == -1) {
	print "\n";
	print "This script extracts sequences with a specific name from all files in a directory and saves them in a separate fasta file.\n";
	print "\nSyntax:\n";
	print " -f X		use file X as input or 'all' to convert all files in the current dir\n";
	print " -n X		search for the name X in sequences from the input file'\n";
	print " -r 		recursively parse all dirs below the current (only usable with '-f all')\n";	
	print "\n";
	exit;
}
my $input = "";
my $recurse = 0;
my $in_ext = "fas";

my $lookup_name = "";

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
	if ($ARGV[$argnum] eq "-n" and  $ARGV[$argnum+1]) {
		$lookup_name = $ARGV[$argnum+1];
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

		my @filelist = $ffor->in(".");
		foreach my $file (@filelist) {
			print "\nExtracting sequences containing $lookup_name from $file...\n";	
			extract_seqs($file,$lookup_name);
		}
	}
	elsif (-e $input) {
		print "\nExtracting sequences containing $lookup_name from $input...\n";	
		extract_seqs($input,$lookup_name);
	}
}
else {
	print "\nNo file specified or file does not exists. Pass a filename after the \"-f\" parameter (or 'all' to use all files in the current dir)\n";
}
exit;

sub extract_seqs {
	my @seqs = ();
	my $file = $_[0];
	my $lookup_name = $_[1];
	my $input = Bio::SeqIO->new(-file => "$file");
	while (my $seq = $input->next_seq){
		my $id = $seq->display_id;
		if ($id =~ $lookup_name) {
			push (@seqs, $seq);
        	}
	}
	my $output = Bio::SeqIO->new(-file => ">$file\_$lookup_name\_extract.fas", -format => "Fasta");
	foreach my $seq_to_write (@seqs) {
		$output->write_seq($seq_to_write);
	}
}

