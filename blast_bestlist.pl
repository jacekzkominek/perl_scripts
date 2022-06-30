#!/usr/bin/perl -w

use strict;
use IO::String;
use Bio::SearchIO;
use File::Find::Object::Rule;

my $input="";

if ($#ARGV == -1) {
	print "\n";
	print "This script prints best hits from BLAST result text files.\n";
	print "\nSyntax:\n";
	print " -f X		use file X as input or 'all' to parse all files in the current dir\n";
	print "\n";
	exit;
}

foreach my $argnum (0 .. $#ARGV) {
	if ($ARGV[$argnum] =~ "-f") {
		$input = $ARGV[$argnum+1];
	}
}
if ($input) {
	if ($input eq "all") {
		my $ffor = File::Find::Object::Rule->file()->name("*.txt");
		$ffor->maxdepth(1);
		my @filelist = $ffor->in(".");
		foreach my $file (@filelist) {
			retrieve_blast_hits($file);
		}
	}
	elsif (-e $input) {
		retrieve_blast_hits($input);
	}
}
else {
	print "\nNo file specified or file does not exists. Pass a filename after the \"-f\" parameter (or 'all' to use all files in the current dir)\n";
}
exit;

sub retrieve_blast_hits {
	my $file = $_[0];
	open (LOG, ">$file\_log");
	print "\nRetrieving hit sequences from $file...\n";

	my $in = new Bio::SearchIO(-format => 'blast', -file   => $file);
	while( my $result = $in->next_result ) {
		while( my $hit = $result->next_hit ) {
			print "\n",$hit->name()."\:\n";
			print LOG "\n",$hit->name()."\:\n";
	    		while( my $hsp = $hit->next_hsp ) {
				if( $hsp->length('total') > 50 ) {
					print "HSP - ",$hsp->score(),"\t",$hsp->start('hit')," - ", $hsp->end('hit'),"\n";
					print LOG "HSP - ",$hsp->score(),"\t",$hsp->start('hit')," - ", $hsp->end('hit'),"\n";
				}
			}
		}
	}
	close (LOG);
	print "\nRetrieval completed.\n";
}
