#!/usr/bin/perl -w

use strict;
use IO::String;
use Bio::AlignIO;
use Bio::SeqIO;
use Bio::SearchIO;
use File::Find::Object::Rule;

my $file_list = "";
my $keyword_list = "";
my $keyword_file = "";

if ($#ARGV == -1) {
	print "\n";
	print "This script sorts results from blast.parse.pl into files with individual hits.\n";
	print "\nSyntax:\n";
	print " -l X		list of words (database names) to look for in the name files\n";
	print " -kw X		list of keywords describing hit names (will be used in the output file)\n";
	print " -kw_f		file with sequences to be used as keywords\n";
	print "\n";
	exit;
}

foreach my $argnum (0 .. $#ARGV) {
	if ($ARGV[$argnum] =~ "-l") {
		$file_list = $ARGV[$argnum+1];
	}
	if ($ARGV[$argnum] =~ "-kw") {
		$keyword_list = $ARGV[$argnum+1];
	}
	if ($ARGV[$argnum] =~ "-kw_f") {
		$keyword_file = $ARGV[$argnum+1];
	}		
}

if ($file_list and ($keyword_file or $keyword_list)) {
	my @files = split("_",$file_list);
	my @kw_list;
	if ($keyword_file ne "") {
		my $seq_in = Bio::SeqIO->new(-file => $keyword_file);
		while (my $seq = $seq_in->next_seq) {
			my $id = $seq->id;
			push (@kw_list, $id);
		}
	}	
	else {
		@kw_list = split("_",$keyword_list);
	}
	foreach my $f (@files) {
		my $ffor = File::Find::Object::Rule->file()->name("*$f*.fas");
		$ffor->maxdepth(1);
		my @filelist = $ffor->in(".");
		foreach my $file (@filelist) {
			print "Parsing file $file...\n";
			foreach my $kw (@kw_list) {
				my $seq_out = Bio::SeqIO->new(-file => ">>$kw.fas", -format => "fasta");
				parse_fasta_file($file,$kw,$seq_out);
			}
		}
	}				
}
exit;

sub parse_fasta_file {
	my $file = $_[0];
	my $kw = $_[1];
	my $seq_out = $_[2];
	my $seq_input = Bio::SeqIO->new(-file => $file);
	while (my $seq = $seq_input->next_seq) {
		if ($seq->id =~ $kw) {
			$seq_out->write_seq($seq);
		}
	}	

}

