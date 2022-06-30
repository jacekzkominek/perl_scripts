#!/usr/bin/perl -w

use strict;
use Bio::Seq;
use Bio::SeqIO;
use Bio::Tools::SeqStats;

my $input = "";
my $seq_obj = Bio::Seq->new();

if ($#ARGV == -1) {
	print "\n";
	print "This script prints molecular weights of protein sequences in a file.\n";
	print "\nSyntax:\n";
	print " -f X		use file X as input\n";
	print "\n";
	exit;
}

foreach my $argnum (0 .. $#ARGV) {
	if ($ARGV[$argnum] eq "-f" and $ARGV[$argnum+1]) {
		$input = $ARGV[$argnum+1];
		if (-e $input) {
			print "\nLoading $input...\n";
		}
		else {
			$input = "";
		}
	}
}

if ($input) {
	if (-e $input) {
		my $seqio_obj = Bio::SeqIO->new(-file => "$input");
		while (my $seq_obj = $seqio_obj->next_seq){   
			my $weight = Bio::Tools::SeqStats->get_mol_wt($seq_obj);
			my $mw="";
			if ($$weight[0] == $$weight[1]) {
				$mw = $$weight[0];
			}
			else {
				$mw = "$$weight[0]"."-"."$$weight[1]";
			}
			print $seq_obj->id(),$seq_obj->desc(),"\t",$mw," Da\n";
		}
	}
}
else {
	print "\nNo file specified or file does not exists. Pass a filename after the \"-f\" parameter.\n";
}
exit;
	
