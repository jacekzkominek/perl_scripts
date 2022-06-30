#!/usr/bin/perl -w

use strict;
use Bio::SeqIO;
use Bio::AlignIO;

my $no_margin="";
my $best_only="";
my $no_extract="";
my $margin = "";
my $ref = "";
my $out = "";
my $margin_size = 300;
my $dblist = "";
my @dbs;
my $sort = 0;

foreach my $argnum (0 .. $#ARGV) {
	if ($ARGV[$argnum] =~ "-f") {
		$ref = $ARGV[$argnum+1];
	}
	if ($ARGV[$argnum] =~ "-out") {
		$out = $ARGV[$argnum+1];
	}
	if ($ARGV[$argnum] =~ "-sort") {
		$sort = 1;
	}
	if ($ARGV[$argnum] =~ "-no_margin") {
		$no_margin = "-no_margin";
	}
	if ($ARGV[$argnum] =~ "-best_only") {
		$best_only = "-best_only";
	}	
	if ($ARGV[$argnum] =~ "-no_extract") {
		$no_extract = "-no_extract";
	}
	if ($ARGV[$argnum] =~ "-margin") {
		$margin = "-margin ";
		$margin_size = $ARGV[$argnum+1];
		$margin = $margin.$margin_size;
	}
	if ($ARGV[$argnum] =~ "-db") {
		if ($ARGV[$argnum+1] =~ "fungi_basic") {
			$dblist = "Acap Afum Anid Aory Ater Bfuc Calb Ccin Cdub Cgla Cglo Cimm Clus Cneo Cpar Ctro Dhan Egos Gmon Gzea Klac Lthe Lklu Lwal Mgri Ncra Pans Pchr Pgui Pnod Rory Sbay Scas Scer Skud Smik Scry Spar Spom Sjap Soct Sscl Tree Umay Uree Vpol Ylip Zrou";
		}
		if ($ARGV[$argnum+1] =~ "fungi_WGD") {
			$dblist = "Cgla Egos Klac Lthe Lklu Lwal Sbay Scas Scer Skud Smik Spar Vpol Zrou";
		}
		if ($ARGV[$argnum+1] =~ "fungi_extended") {
			$dblist = "Acap Afum Anid Aory Ater Bfuc Calb Ccin Cdub Cgla Cglo Cimm Clus Cneo Cpar Ctro Dhan Egos Gmon Gzea Klac Lthe Lklu Lwal Mgri Ncra Pans Pchr Pgui Pnod Rory Sbay Scas Scer Skud Smik Scry Spar Spom Sjap Soct Sscl Tree Umay Uree Vpol Ylip Zrou";
		}
		if ($ARGV[$argnum+1] =~ "4clades") {
			$dblist = "Acap Afum Anid Aory Ater Bfuc Calb Cdub Cgla Cglo Cimm Clus Cpar Ctro Dhan Egos Gmon Gzea Klac Lklu Lwal Mgri Ncra Pans Pgui Pnod Sbay Scas Scer Skud Smik Smik_WashU Spar Sscl Tree Uree";
		}
		if ($ARGV[$argnum+1] =~ "archea") {
			$dblist = "Aboo Aful Aper Apro Asac Cmaq Dkam Fpla Hbor Hbut Hlac Hmar Hmuk Hsal Htur Huta Hvol Hwal Iagg Ihos Kcry Maeo Mbar Mboo Mbur Meve Mfer Mhun Minf Mjan Mkan Mlab Mmah Mmar Mmaz Moki Mpal Mpet Mrum Msed Msmi Msp. Msta Mthe Mvan Mvol Mvul Nequ Nmag Nmar Npha Paby Paer Pars Pcal Pfur Phor Pisl Pkod Ptor Saci Shel Sisl Smar Ssol Stok Taci Tagg Tgam Tneu Tonn Tpen Tsib Tvol Umet Vdis";
		}
		if ($ARGV[$argnum+1] =~ "basidio") {
			$dblist = "Adel Cput Dspe Dsqu Fmed Gluc Lbic Mglo Mlar Mvio Mosm Pgra Ptri Pstr Scom Slac Shir Tver Tmes Wseb";
		}
		if ($ARGV[$argnum+1] =~ "other") {
			$dblist = "Amac Bden Spun";
		}				
		@dbs = split (" ",$dblist);		
		if ($ARGV[$argnum+1] =~ "_aa") {
			foreach my $db (@dbs) {
				$db = $db."_aa";
			}
		}
		elsif ($ARGV[$argnum+1] =~ "_nn") {
			foreach my $db (@dbs) {
				$db = $db."_nn";
			}
		}
	}
}

foreach my $db (@dbs) {
	system "blast_parse.pl -f $ref -db $db -out $out $no_margin $best_only $no_extract $margin";
}
if ($sort == 1) {
#	my $seq_in = Bio::SeqIO->new(-file => $ref);
#	my $kw_list="";
#	$dblist =~ s/\s/_/g; 
#	while (my $seq = $seq_in->next_seq) {
#		my $id = $seq->id;
#		$kw_list = $kw_list.$id.("_");
#	}
#	$kw_list = substr($kw_list, 0, -1);
#
#	print "blast_sort.pl -l $dblist -kw $kw_list";
	
	print "blast_sort.pl -l $dblist -kw_f $ref";

}


exit;

