#!/usr/bin/perl -w

use strict;
use Bio::Seq;
use Bio::SeqIO;
use Cwd;
use File::Find::Object::Rule;

my $dir = getcwd();
my $out = "horiz_preds.fas";
my $count = 0;

if ($#ARGV == -1) {
	print "\n";
	print "This script joins Psipred secondary structure predictions (*.horiz files) from the selected directory into a single fasta file.\n";
	print "\nSyntax:\n";
	print " -dir X		load files from directory X\n";
	print " -out X		output to file X\n";
	print "\n";
	exit;
}

foreach my $argnum (0 .. $#ARGV) {
	if ($ARGV[$argnum] eq "-dir") {
		$dir = $ARGV[$argnum+1];
	}
	if ($ARGV[$argnum] eq "-out") {
		$out = $ARGV[$argnum+1];
	}
}
my $out_seqfile = Bio::SeqIO->new(-file => ">$out");

print "\nFinding Psipred *.horiz files in \"$dir\" and writing them out to $out...";
my $ffor = File::Find::Object::Rule->file()->name("*.horiz");
$ffor->maxdepth(1);
my @filelist = $ffor->in($dir);
foreach my $file (@filelist) {
	extract_psipred_pred($file,$out_seqfile);
}
print "done. Located $count files.\n\n";

exit;

sub extract_psipred_pred {
	my $file =  $_[0];
	my $outfile = $_[1];
	my $pred = "";
	open (HORIZ_FILE, $file);
	while (<HORIZ_FILE>) {
		chomp;
		my @fields = split(' ');	
		my $fields_count = @fields;
		if ($fields_count >= 1 and $fields[0] eq "Pred:") {
			$pred = $pred.$fields[1];
		}
	}
	my $seq = Bio::Seq->new( -seq => $pred, -id  => $file);
	$outfile->write_seq($seq);
	$count++;
}




