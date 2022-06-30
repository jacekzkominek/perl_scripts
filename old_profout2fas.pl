#!/usr/bin/perl -w

use strict;
use Bio::Seq;
use Bio::SeqIO;
use Cwd;
use File::Find::Object::Rule;

my $dir = getcwd();
my $out = "prof_preds.fas";
my $count = 0;

if ($#ARGV == -1) {
	print "\n";
	print "This script joins prof secondary structure predictions (*.out files) from the selected directory into a single fasta file.\n";
	print "\nSyntax:\n";
	print " -dir X		load files from directory X (use \"\.\" for current directory)\n";
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

print "\nFinding Prof *.out files in \"$dir\" and writing them out to $out...";
my $ffor = File::Find::Object::Rule->file()->name("*.out");
$ffor->maxdepth(1);
my @filelist = $ffor->in($dir);
foreach my $file (@filelist) {
	extract_prof_pred($file,$out_seqfile);
}
print "done. Located $count files.\n\n";

exit;


sub extract_prof_pred {
	my $file =  $_[0];
	my $outfile = $_[1];
	my $pred = "";

	open (PROFOUT_FILE, $file);
	while (<PROFOUT_FILE>) {
		chomp;
		my @fields = split(' ');
		my $fields_count = 0; 
		$fields_count = @fields;
		if ($fields_count >= 21 and ($fields[21] eq "C" or $fields[21] eq "E" or $fields[21] eq "H")) {
			$pred = $pred.$fields[21];
		}
	}
	my $seq = Bio::Seq->new( -seq => $pred, -id  => $file);
	$outfile->write_seq($seq);
	$count++;
}

