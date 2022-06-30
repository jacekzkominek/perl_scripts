#!/usr/bin/perl

use strict;
use File::Copy;
use File::Find::Object::Rule;

my $input = "";
my $prot = "F";

if ($#ARGV == -1) {
	print "\n";
	print "This script creates BLAST databases for all Fasta (*.fas) files in the current directory. It expects the file names to be in the form: \"genus-name_species-name.fas\" and will truncate the resulting database name to the first letter of the genus and first 3 letters of the species name.\n";
	print "\nSyntax:\n";
	print " -f X		use file X as input or 'all' to parse all files in the current dir\n";
	print " -prot		use this switch if protein sequences are used\n";
	print "\n";
	exit;
}

foreach my $argnum (0 .. $#ARGV) {
	if ($ARGV[$argnum] =~ "-f") {
		$input = $ARGV[$argnum+1];
	}
	if ($ARGV[$argnum] =~ "-prot") {
		$prot = "T";
	}	
}
if ($input) {
	if ($input eq "all") {
		my $ffor = File::Find::Object::Rule->file()->name("*.fas");
		$ffor->maxdepth(1);
		my @filelist = $ffor->in(".");
		foreach my $file (@filelist) {
			format_seq_files($file);
		}
	}
	elsif (-e $input) {
		format_seq_files($input);
	}
}
else {
	print "\nNo file specified or file does not exists. Pass a filename after the \"-f\" parameter (or 'all' to use all files in the current dir)\n";
}
exit();

sub format_seq_files {
	my $file = $_[0];
	print "$file\n";
	my $spec1 = substr($file,0,1);
	my $spec2 = substr($file,index($file,"\_")+1,3);
	my $spec = ucfirst($spec1.$spec2);
	my $dbtype = "nucl";
	if ($prot eq "T") {
		$spec = $spec."\_aa";
		$dbtype = "prot";
	}
	else {
		$spec = $spec."\_nn";		
	}
	my $run_dbf = "makeblastdb -in $file -title $spec -dbtype $dbtype -out $spec";
	system $run_dbf;
	copy($file, "$spec\.fas");
}

