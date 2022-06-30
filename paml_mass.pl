#!/usr/bin/perl -W

use strict;
use Bio::Tools::Run::Phylo::PAML::Codeml;
use Bio::AlignIO;
use Cwd;

if ($#ARGV == -1) {
	print "\n";
	print "This script runs PAML's codeml on a sequence file+tree file in a pipeline with variety of options.\n";
	print "\nSyntax:\n";
	print " -f X		use file X as input\n";
	print " -treef X	use treefile X as input\n";
	print " -out X		output results to file X\n";
	print " -model X_Y...	Use different models X,Y etc. (separated by the underscore)\n";
	print " -NSsites X_Y...	use different NSsites values X,Y etc. (separated with the underscore) \n";
	print " -ncatG X_Y...	use different amounts of gamma categories (separated with the underscore)\n";	
	print " -codonfreq X	use CodonFreq value X\n";
	print " -kappa X_Y...	estimate kappa from different starting values X,Y etc. (separated with the underscore)\n";
	print " -kappaf X_Y...	use different fixed kappa values X,Y etc. (separated with the underscore)\n";
	print " -omega X_Y...	estimate omega from different starting values X,Y etc. (separated with the underscore)\n";
	print " -omegaf X_Y...	use different fixed omega values X,Y etc. (separated with the underscore)\n";
	print " -no_getSE	do not get standard errors of estimates\n";
	print " -no_clean	do not remove ambiguity characters\n";
	print " -verbose	print codeml output on the screen once a run is finished\n";
	print " -log		print codeml output to a log file once a run is finished\n";	
	print " -fix_blength X	Use fix_blength value X (1 by default)\n";		
	print "\n";
	exit;
}

my $input = "";
my $output = "";
my $tree_input = "";
my $model = "0_";
my $NSsites = "0_";
my $codfreq = 3;
my $ncatG = "1_";
my $kappa = "2_";
my $kappaf = 0;
my $omega = "1_";
my $omegaf = 0;
my $no_getSE = 0;
my $no_clean = 0;
my $verbose = 0;
my $log = 0;
my $dnds_sites = 0;
my $fix_blength = 1;

foreach my $argnum (0 .. $#ARGV) {
	if ($ARGV[$argnum] eq "-f" and $ARGV[$argnum+1]) {
		$input = $ARGV[$argnum+1];
		if (-e $input) {
			print "\nLoading $input...\n";
		}
		else {
			$input="";
		}
	}
	if ($ARGV[$argnum] eq "-out") {
		$output = $ARGV[$argnum+1];
	}
	if ($ARGV[$argnum] eq "-treef" and -e $ARGV[$argnum+1]) {
		$tree_input = $ARGV[$argnum+1];
	}
	if ($ARGV[$argnum] eq "-model") {
		$model = $ARGV[$argnum+1];
	}	
	if ($ARGV[$argnum] eq "-NSsites") {
		$NSsites = $ARGV[$argnum+1];
	}
	if ($ARGV[$argnum] eq "-ncatG") {
		$ncatG = $ARGV[$argnum+1];
	}
	if ($ARGV[$argnum] eq "-codfreq") {
		$codfreq = $ARGV[$argnum+1];
	}	
	if ($ARGV[$argnum] eq "-kappa") {
		$kappa = $ARGV[$argnum+1];
	}
	if ($ARGV[$argnum] eq "-kappaf") {
		$kappaf = 1;
	}	
	if ($ARGV[$argnum] eq "-omega") {
		$omega = $ARGV[$argnum+1];
	}					
	if ($ARGV[$argnum] eq "-omegaf") {
		$omegaf = 1;
	}
	if ($ARGV[$argnum] eq "-no_getSE") {
		$no_getSE = 1;
	}
	if ($ARGV[$argnum] eq "-no_clean") {
		$no_clean = 1;
	}	
	if ($ARGV[$argnum] eq "-verbose") {
		$verbose = 1;
	}							
	if ($ARGV[$argnum] eq "-log") {
		$log = 1;
	}
	if ($ARGV[$argnum] eq "-dnds_sites") {
		$dnds_sites = 1;
	}		
	if ($ARGV[$argnum] eq "-fix_blength") {
		$fix_blength = $ARGV[$argnum+1];
	}		
}	
}

my $alignio = Bio::AlignIO->new(-format => 'fasta', -file => $input);
my $aln = $alignio->next_aln;
my $treeio = Bio::TreeIO->new(-format => 'newick', -file => $tree_input);
my $tree = $treeio->next_tree;
my $false = 0;
my $codeml = Bio::Tools::Run::Phylo::PAML::Codeml->new(-alignment => $aln, -tree => $tree, -save_tempfiles => 0);

#SET FIXED PARAMETERS AND DEFAULTS
$codeml->set_parameter('noisy',3);
$codeml->set_parameter('runmode',0);
$codeml->set_parameter('CodonFreq',$codfreq);
$codeml->set_parameter('clock',0);
$codeml->set_parameter('seqtype',1);
$codeml->set_parameter('icode',0);
$codeml->set_parameter('Mgene',0);
$codeml->set_parameter('cleandata',1-$no_clean);
$codeml->set_parameter('RateAncestor',0);
$codeml->set_parameter('fix_blength',$fix_blength);
$codeml->set_parameter('getSE',1-$no_getSE);
$codeml->set_parameter('model',0);
$codeml->set_parameter('NSsites',0);
$codeml->set_parameter('ncatG',1);
$codeml->set_parameter('fix_kappa',$kappaf);
$codeml->set_parameter('kappa',2);
$codeml->set_parameter('fix_omega',$omegaf);
$codeml->set_parameter('omega',1);


#RUN WITH CUSTOM PARAMETERS
my @ms = split('_', $model);
my @ns = split('_', $NSsites);
my @cats = split('_', $ncatG);
my @ks = split('_', $kappa);
my @os = split('_', $omega);

$codeml->no_param_checks(1);
foreach my $m (@ms) {
	$codeml->set_parameter('model',$m);
	foreach my $n (@ns) {
		$codeml->set_parameter('NSsites',$n);
		foreach my $cat (@cats) {
			$codeml->set_parameter('ncatG',$cat);
			foreach my $k (@ks) {
				$codeml->set_parameter('kappa',$k);
				foreach my $o (@os) {
					$codeml->set_parameter('omega',$o);
					my $outfile_name = getcwd."\/".$output."M$m"."_N$n"."_Cat$cat"."_cf$codfreq"."_K$k"."_O$o";
					$codeml->outfile_name($outfile_name);
					print "\nRunning codeml on $input with model=$m, NS=$n, ncatG=$cat, codonFreq=$codfreq, kappa=$k, omega=$o\n";
					my ($rc,$parser) = $codeml->run();
					if ($dnds_sites == 1) {
						my $result = $parser->next_result;
						my $classes = $result->num_site_classes();
						if ($classes > 1) {
							my $dNdS = $result->dnds_site_classes();
							print "\n---\n";	
							for (my $i1 = 1 ; $i1 <= $classes ; $i1++) {
								print $dNdS->{w}."\n";
							}
							print "\n---\n";					
						}
					}
					if ($verbose == 1) {
						print("codeml error : ", $codeml->error_string, "\n");
					}
					if ($log == 1) {
						open (LOG, ">$outfile_name\_log");
						print LOG ($codeml->error_string,"\n");
						close (LOG);
					}
				}
			}
		}		
		
	}		
}


exit();

