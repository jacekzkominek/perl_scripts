#!/usr/bin/perl -w

use strict;
use Bio::Seq;
use Bio::SeqIO;
use Cwd;

my $input = "";

if ($#ARGV == -1) {
	print "\n";
	print "This script parses and summarizes the output from the SignalP program for multiple sequences.\n";
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
	my $seq_name = "";
	my $name = "";
	my $count = 0;
	my $min_nn_pos = 9999;
	my $min_hmm_pos = 9999;
	my $max_nn_pos = 0;
	my $max_hmm_pos = 0;
	open(LOG, ">$input\_signalp_log.txt");
	open(SIGNALP_FILE, $input);
	while (<SIGNALP_FILE>) {
		chomp;
		my @fields = split(' ');	
		if ($fields[0]) {
			if ($fields[0] =~ ">") {
				$name = substr($fields[0],1,length($fields[0]));
				if ($name ne $seq_name) {
					$seq_name = $name;
					print "\n$seq_name\n";
					print LOG "\n$seq_name\n";
					$count++;
				}
			}
			if ($fields[0] eq "#" and $fields[1] eq "Most") {
				print "NN prediction: $fields[7]-$fields[9] ($fields[10])\n";		
				print LOG "NN prediction: $fields[7]-$fields[9] ($fields[10])\n";		
				if ($fields[7] < $min_nn_pos) {
					$min_nn_pos = $fields[7];
				}
				if ($fields[7] > $max_nn_pos) {
					$max_nn_pos = $fields[7];
				}
			}
			if ($fields[0] eq "Max") {
				print "HMM prediction ($fields[4]): $fields[7]-$fields[9]\n";		
				print LOG "HMM prediction ($fields[4]): $fields[7]-$fields[9]\n";
				if ($fields[7] < $min_hmm_pos) {
					$min_hmm_pos = $fields[7];
				}
				if ($fields[7] > $max_hmm_pos) {
					$max_hmm_pos = $fields[7];
				}
			}
		}
	}
	print "\n";
	print LOG "\n";

	print "NN cleavage position range: $min_nn_pos-$max_nn_pos\n";
	print LOG "NN cleavage position range: $min_nn_pos-$max_nn_pos\n";
	print "HMM cleavage position range: $min_hmm_pos-$max_hmm_pos\n";
	print LOG "HMM cleavage position range: $min_hmm_pos-$max_hmm_pos\n";

	close (SIGNALP_FILE);
	close (LOG);
}
else {
	print "\nNo file specified or file does not exists. Pass a filename after the \"-f\" parameter.\n";
}

exit;
