#!/usr/bin/perl -w

use strict;
use File::Find::Object::Rule;

my $input = "";
my $output = "codeml.log";
my $verbose = 0;
if ($#ARGV == -1) {
	print "\n";
	print "This script parses output from the codeml program and extracts the LnL, dN/dS and kappa values .\n";
	print "\nSyntax:\n";
	print " -f X		use file X as input or 'all' to convert all files in the current dir\n";
	print " -out X		output to file X\n";
	print " -verbose	show output for all parsed files";
	print "\n";
	exit;
}

foreach my $argnum (0 .. $#ARGV) {
	if ($ARGV[$argnum] eq "-f" and $ARGV[$argnum+1]) {
		$input = $ARGV[$argnum+1];
	}
	if ($ARGV[$argnum] eq "-out" and $ARGV[$argnum+1]) {
		$output = $ARGV[$argnum+1];
	}
	if ($ARGV[$argnum] eq "-verbose") {
		$verbose = 1;
	}	
}

my @lnls = ();
my @kappas = ();
my @lnls_kappas = ();
my $file_count=0;
if ($input) {
	open (LOG, ">$output");
	if ($input eq "all") {
		my $ffor = File::Find::Object::Rule->file()->name("*");
		$ffor->maxdepth(1);

		my @filelist = $ffor->in(".");
		foreach my $file (@filelist) {
			$file_count++;
			parse_paml_file($file);
		}
	}
	elsif (-e $input) {
		parse_paml_file($input);
	}
}
else {
	print "\nNo file specified. Use the \"-f\" param to pass a filename or 'all'.\n";
}
print "\nInfo from all analyzed runs, sorted by log-likelihood:";
print LOG "\nInfo from all analyzed runs, sorted by log-likelihood:";

@lnls_kappas = sort @lnls_kappas;
print "\n",@lnls_kappas,"\n";
print LOG "\n",@lnls_kappas,"\n";

@lnls = sort @lnls;
@kappas = sort @kappas;
#print "\n",@lnls,"\n";
#print "\n",@kappas,"\n";

print LOG "\n",@lnls,"\n";
print LOG "\n",@kappas,"\n";

close LOG;
exit;

sub parse_paml_file {
	my $input = $_[0];
	open (PAML_FILE, $input);
	
	if ($verbose == 1) {
		print "\nLoading file $input...\n";
		print LOG "\nLoading file $input...\n";
	}
	
	my $lnl = 0;
	my $kappa = 0;
	my $omega_block = 0;
	my $classes = 0;
	my @data = ();
	while (<PAML_FILE>) {
		chomp;
		if ($_ eq "") {
			next;
		}
		if ($_ =~ "lnL") {
			my @line = split(' ');
			foreach my $l (@line) {
				if (substr($l,0,1) eq '-') {
					$lnl = $l;
					my $lnl_2f = sprintf("%.2f", $lnl);
					push (@lnls, $lnl_2f."_".$input."\n");		
				}
			}
			next;
		}
		if ($_ =~ "kappa") {
			my @line = split(' ');
			$kappa = $line[3];
			my $kappa_2f = sprintf("%.2f", $kappa);
			push (@kappas, $kappa_2f."_".$input."\n");
			push (@lnls_kappas, sprintf("%.2f", $lnl)."___K".$kappa_2f."___".$input."\n");
			next;
		}
		if ($_ =~ "dN/dS " and $_ =~ "site classes") {
			$omega_block = 1;
			my @line = split(' ');
			$classes = $line[5];
			$classes = substr($classes,3,1);
			next;
		}
		if ($omega_block == 1) {
			if ($_ =~ "Bayes Empirical Bayes") {
				$omega_block = 0;
				next;
			}
			if ($classes > 0) {
				my @line = split(' ');
				push (@data, grep (/[0-9]+\.[0-9]+/, @line));
				next;
			}
		}
	}
	if ($classes > 0 and $verbose == 1) {
		my $data_size = @data;
		my $types = ($data_size/$classes)-1;
		print "\n";
		print "LnL=",$lnl,"\n";
		print "Kappa=",$kappa,"\n";
		print "Site classes=",$classes,"\n";
		print "Proportion:\t";

		print LOG "\n";
		print LOG "LnL=",$lnl,"\n";		
		print LOG "Kappa=",$kappa,"\n";		
		print LOG "Site classes=",$classes,"\n";
		print LOG "Proportion:\t";

		for (my $i1 = 0; $i1 < $classes ;$i1++) {
			print $data[$i1]."\t";
			print LOG $data[$i1]."\t";		
		}
		for (my $i1 = $classes, my $br = 1; $i1 < $data_size ;$i1++) {
			if ($i1 % $classes == 0) {
				print "\nBranch $br:\t";
				print LOG "\nBranch $br:\t";			
				$br++;
			}
			print $data[$i1]."\t";
			print LOG $data[$i1]."\t";		
		}
		print "\n";
		print LOG "\n";
	}
	close PAML_FILE;
}


