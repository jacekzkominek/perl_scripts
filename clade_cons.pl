#!/usr/bin/perl -w

use strict;
use IO::String;
use Bio::Seq; 
use Bio::SeqIO;
use Bio::AlignIO;
use File::Find::Object::Rule;


my $input = "";
my $clade_split = 1;
my $align_split_pos = 0;
my $seq_split_pos = 0;
my $seq_split_ref = 1;
my $strict_clade1 = 0;
my $strict_clade2 = 0;
my $only_count = 0;
my $print_counts = 0;
my $tolerance = 0;

if ($#ARGV == -1) {
	print "\n";
	print "This script identifies conserved positions in an alignment. Positions conserved with clade-specific conservation between 2 groups of sequences.\n";
	print "\nSyntax:\n";
	print " -f X		use file X as input\n";
	print " -clade_split X	use all sequences before X-th as the first group, and the rest (including X-th) as the second group.\n";
	print " -align_split X  separately calculate conserved positions before and after position X in the alignment.\n";
	print " -seq_split X_Y	separately calculate conserved positions before and after position X in sequence Y (if Y is skipped, first sequence is taken.\n";
	print " -strict_clade1	print positions conserved only in the first clade.\n";
	print " -strict_clade2	print positions conserved only in the second clade.\n";
	print " -print_counts	print the quantities of non-conserved residues as total counts instead of proportions.\n";
	print " -only_count	print only a summary, without individual conserved positions.\n";
	print " -tolerance X	allow X non-conforming residues in the annotation process.\n";
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
	if ($ARGV[$argnum] eq "-clade_split") {
		$clade_split = $ARGV[$argnum+1];
	}
	if ($ARGV[$argnum] eq "-align_split") {
		$align_split_pos = $ARGV[$argnum+1];
	}	
	if ($ARGV[$argnum] eq "-seq_split") {
		my $arg = $ARGV[$argnum+1];
		$seq_split_pos = substr($arg, 0, rindex($arg, "_"));
		$seq_split_ref = substr($arg, rindex($arg, "_")+1);
	}	
	if ($ARGV[$argnum] eq "-strict_clade1") {
		$strict_clade1 = 1;
	}
	if ($ARGV[$argnum] eq "-strict_clade2") {
		$strict_clade2 = 1;
	}
	if ($ARGV[$argnum] eq "-only_count") {
		$only_count = 1;
	}
	if ($ARGV[$argnum] eq "-print_counts") {
		$print_counts = 1;
	}							
	if ($ARGV[$argnum] eq "-tolerance") {
		$tolerance = $ARGV[$argnum+1];
	}	
}

if ($input) {
	open (LOG, ">$input\_log");

	print "\nLoading sequences from $input...\n";
	print LOG "\nLoading sequences from $input...\n";
	
	my $seq_in = Bio::AlignIO->new(-format => 'fasta', -file => $input);
	my $seq_aln = $seq_in->next_aln;
	my $seq_count = $seq_aln->num_sequences;
	my $aln_len = $seq_aln->length;
	if ($seq_count < 2) {
		print "\nLess than 2 sequences in the file. Exiting\n";
		exit;
	}	
	
	my $clade1 = new Bio::SimpleAlign();
	my $clade2 = new Bio::SimpleAlign();
	my @univ_cons_res_arr1 = ();
	my @univ_cons_res_arr2 = ();
	my @ref_cons_res_arr1 = ();
	my @ref_cons_res_arr2 = ();
	my @univ_cons_arr = ();
	my @univ_cons_short1 = ();
	my @univ_cons_short2 = ();

	if ($clade_split == 1) {
		#Calculate universally conserved positions
		my $seq1 = $seq_aln->get_seq_by_pos(1);	
		my $seq2 = $seq_aln->get_seq_by_pos(2);
		my @temp_arr = @{get_univ_cons($seq1, $seq2)};
		if ($seq_count == 2) {
			@univ_cons_arr = @temp_arr;
		}
		my @univ_cons_arr = ();
		foreach my $pos (@temp_arr) {
			for (my $i1 = 2; $i1 <= $seq_count; $i1++) {
				my $seq2 = $seq_aln->get_seq_by_pos($i1);
				my $res1 = $seq1->subseq($pos,$pos);
				my $res2 = $seq2->subseq($pos,$pos);
				if ($res1 ne $res2) {
					last;
				}
				if ($i1 == $seq_count) {
					push (@univ_cons_arr, $pos);
				}
			}
		}
	}
	else {
		#Calculate clade-specifically conserved positions
		if ($seq_count <= 2) {
			print "\nOnly 2 sequences in the file, exiting.\n\n";
			print LOG "\nOnly 2 sequences in the file, exiting.\n\n";
			exit;
		}
	}
			
	for (my $i1 = 1; $i1 <= $seq_count; $i1++) {
		my $seq = $seq_aln->get_seq_by_pos($i1);
		if ($i1 < $clade_split) {
			$clade1->add_seq($seq);
		}
		else {
			$clade2->add_seq($seq);
		}
	}
	my $clade1_count = $clade1->num_sequences;
	my $clade2_count = $clade2->num_sequences;
	
	
	print "First clade has $clade1_count sequences:\n";
	print LOG "First clade has $clade1_count sequences:\n";
	foreach my $seq1 ($clade1->each_seq()) {
		print $seq1->id(),"\n";
		print LOG $seq1->id(),"\n";
	} 
	print "\nSecond clade has $clade2_count sequences:\n";
	print LOG "\nSecond clade has $clade2_count sequences:\n";
	foreach my $seq2 ($clade2->each_seq()) {
		print $seq2->id(),"\n";
		print LOG $seq2->id(),"\n";
	}
	print "\nFinding conserved positions within clades...\n";
	print LOG "\nFinding conserved positions within clades...\n";

	#PAIRWISE CHECK OF CLADE 1	
	my $ref_seq1 = $clade1->get_seq_by_pos(1);	
	for (my $i1 = 2; $i1 <= $clade1_count; $i1++) {
		my $seq1 = $clade1->get_seq_by_pos($i1);
		my @temp_arr = @{get_univ_cons($ref_seq1, $seq1)};
		if ($i1 == 2) {
			push (@ref_cons_res_arr1, @temp_arr);		
		}
		push (@univ_cons_res_arr1, @temp_arr);
	}
	
	#PAIRWISE CHECK OF CLADE 2
	my $ref_seq2 = $clade2->get_seq_by_pos(1);
	#print $ref_seq2->id();
	for (my $i2 = 2; $i2 <= $clade2_count; $i2++) {
		my $seq2 = $clade2->get_seq_by_pos($i2);
		my @temp_arr = @{get_univ_cons($ref_seq2, $seq2)};
		if ($i2 == 2) {
			push (@ref_cons_res_arr2, @temp_arr);		
		}
		push (@univ_cons_res_arr2, @temp_arr);
	}

	#FINAL ARRAY OF CONSERVED POSITIONS OF CLADE1
	foreach my $ref_res1 (@ref_cons_res_arr1) {
		my $res_count1 = 1;
		foreach my $seq_res1 (@univ_cons_res_arr1) {
			if ($ref_res1 == $seq_res1 or $ref_res1 eq "X" or $seq_res1 eq "X") {
				$res_count1++;
			}
		}
		if ($res_count1 == $clade1_count) {
			push (@univ_cons_short1, $ref_res1);
		}
	}

	#FINAL ARRAY OF CONSERVED POSITIONS OF CLADE2
	foreach my $ref_res2 (@ref_cons_res_arr2) {
		my $res_count2 = 1;
		foreach my $seq_res2 (@univ_cons_res_arr2) {
			if ($ref_res2 == $seq_res2 or $ref_res2 eq "X" or $seq_res2 eq "X") {
				$res_count2++;
			}
		}
		if ($res_count2 == $clade2_count) {
			push (@univ_cons_short2, $ref_res2);
		}
	}	

	my $c1_size = @univ_cons_short1;
	my $c2_size = @univ_cons_short2;
	print "Positions conserved universally in clade 1: $c1_size \/ $aln_len \n";
	print LOG "Positions conserved universally in clade 1: $c1_size \/ $aln_len \n";
	print "Positions conserved universally in clade 2: $c2_size \/ $aln_len \n";
	print LOG "Positions conserved universally in clade 2: $c2_size \/ $aln_len \n";

	#CHECK FOR POSITIONS CONSERVED IN BOTH CLADES - QUICK METHOD
#	my @array1 = @univ_cons_short1;
#	my @array2 = @univ_cons_short2;
#	my %in_array2 = map { $_ => 1 } @array2;
#	my @array3 = grep { $in_array2{$_} } @array1;
#	print "common @array3\n";
#	my $int_size = @array3;
#	print "Positions identical in both clades: $int_size \/ $aln_len \n";
#	print LOG "Positions identical in both clades: $int_size \/ $aln_len \n";
	
	#CHECK FOR POSITIONS CONSERVED IN BOTH CLADES - FULL METHOD
	if ($only_count == 0) {
		print "\nPositions conserved in both clades:\n";
		print LOG "\nPositions conserved in both clades:\n";
	}
	my $univ_cons_count = 0;
	for (my $i1=1; $i1 < $c1_size; $i1++) {
		my $c1_cons_pos = $univ_cons_short1[$i1];
		for (my $i2=1; $i2 < $c2_size; $i2++) {
			my $c2_cons_pos = $univ_cons_short2[$i2];
			if ($c1_cons_pos == $c2_cons_pos) {
				my $res1 = $ref_seq1->subseq($c1_cons_pos,$c1_cons_pos);
				my $res2 = $ref_seq2->subseq($c2_cons_pos,$c2_cons_pos);
				if ($res1 eq $res2) {
					if ($only_count == 0) {					
						print "Position $c1_cons_pos: $res1\n"; 
						print LOG "Position $c1_cons_pos: $res1\n"; 
					}
					$univ_cons_count++;
					last;
				}
			}
		}
	}	
	print "Positions conserved in both clades: $univ_cons_count \/ $aln_len \n";
	print LOG "Positions conserved in both clades: $univ_cons_count \/ $aln_len \n";
		
	if ($only_count == 0) {	
		print "\nPositions conserved clade-specifically (type-II sites) and their EI rank:\n";
		print LOG "\nPositions conserved clade-specifically (type-II sites) and their EI rank:\n";
	}
	my $clade_spec_cons_count = 0;
	for (my $i1=1; $i1 < $c1_size; $i1++) {
		my $c1_cons_pos = $univ_cons_short1[$i1];
		for (my $i2=1; $i2 < $c2_size; $i2++) {
			my $c2_cons_pos = $univ_cons_short2[$i2];
			if ($c1_cons_pos == $c2_cons_pos) {
				my $res1 = $ref_seq1->subseq($c1_cons_pos,$c1_cons_pos);
				my $res2 = $ref_seq2->subseq($c2_cons_pos,$c2_cons_pos);
				if ($res1 ne $res2) {					
					my $ei_pos = "";
					if (check_EI_pos($res1.$res2) > 0) {
						$ei_pos = "#".check_EI_pos($res1.$res2);
					}
					else {
						$ei_pos = "multi";
					}
					if ($only_count == 0) {					
						print "Position $c1_cons_pos:\t$res1->$res2\t($ei_pos)\n"; 
						print LOG "Position $c1_cons_pos:\t$res1->$res2\t($ei_pos)\n"; 
					}
					$clade_spec_cons_count++;
					last;
				}
			}
		}
	}
	print "Positions conserved clade-specifically (type-II sites): $clade_spec_cons_count\n";
	print LOG "Positions conserved clade-specifically (type-II sites): $clade_spec_cons_count\n";



	if ($strict_clade1 == 1) {
		if ($only_count == 0) {
			print "\nPositions conserved only in the first clade (type-I sites):\n";	
			print LOG "\nPositions conserved only in the first clade (type-I sites):\n";	
		}
		my $c1_spec_cons_count = 0;
		for (my $i1=1; $i1 < $c1_size; $i1++) {
			my $c1_cons_pos = $univ_cons_short1[$i1];
			for (my $i2=1; $i2 <= $c2_size; $i2++) {
				if ($i2 == $c2_size) {
					my $res1 = $ref_seq1->subseq($c1_cons_pos,$c1_cons_pos);
					my $ref_seq2 = $clade2->get_seq_by_pos(1);
					my %res_count = ();
					$res_count{$ref_seq2->subseq($c1_cons_pos,$c1_cons_pos)} = 1;
					for (my $i3 = 2; $i3 <= $clade2_count; $i3++) {
						my $seq2 = $clade2->get_seq_by_pos($i3);
						my $seq2_res = $seq2->subseq($c1_cons_pos,$c1_cons_pos);
						if (!(defined $res_count{$seq2_res})) {
							$res_count{$seq2_res} = 1;
						}
						elsif ($res_count{$seq2_res} > 0) {
							$res_count{$seq2_res} = $res_count{$seq2_res}+1;
						}
					}
					my $c2_residues = "";
					foreach my $k (keys (%res_count)) {
						if ($print_counts == 0) {
							my $freq = sprintf("%.3f", $res_count{$k}/$clade2_count);
							$c2_residues = $c2_residues.$k."-".($freq*100)."% ";
						}
						else {
							$c2_residues = $c2_residues.$k."-".$res_count{$k}." ";
						}
					}
					$c2_residues =~ s/\s+$//;
					if ($only_count == 0) {
						print "Position $c1_cons_pos: $res1 (clade2: $c2_residues)\n"; 
						print LOG "Position $c1_cons_pos: $res1 (clade2: $c2_residues)\n"; 
					}
					$c1_spec_cons_count++;
					last;			
				}
				my $c2_cons_pos = $univ_cons_short2[$i2];
				if ($c1_cons_pos == $c2_cons_pos) {
					last;
				}
			}
		}	
		print "Positions conserved in the first clade (type-I sites): $c1_spec_cons_count\n";
		print LOG "Positions conserved in the first clade (type-I sites): $c1_spec_cons_count\n";
	} 

	if ($strict_clade2 == 1) {	
		if ($only_count == 0) {
			print "\nPositions conserved only in the second clade (type-I sites):\n";
			print LOG "\nPositions conserved only in the second clade (type-I sites):\n";
		}
		my $c2_spec_cons_count = 0;		
		for (my $i1=1; $i1 < $c2_size; $i1++) {
			my $c2_cons_pos = $univ_cons_short2[$i1];
			for (my $i2=1; $i2 <= $c1_size; $i2++) {
				if ($i2 == $c1_size) {
					my $res2 = $ref_seq2->subseq($c2_cons_pos,$c2_cons_pos);
					my $ref_seq1 = $clade1->get_seq_by_pos(1);

					my %res_count = ();
					$res_count{$ref_seq1->subseq($c2_cons_pos,$c2_cons_pos)} = 1;
					for (my $i3 = 2; $i3 <= $clade1_count; $i3++) {
						my $seq2 = $clade1->get_seq_by_pos($i3);
						my $seq2_res = $seq2->subseq($c2_cons_pos,$c2_cons_pos);
						if (!(defined $res_count{$seq2_res})) {
							$res_count{$seq2_res} = 1;
						}
						elsif ($res_count{$seq2_res} > 0) {
							$res_count{$seq2_res} = $res_count{$seq2_res}+1;
						}
					}
					my $c1_residues = "";
					foreach my $k (keys (%res_count)) {
						if ($print_counts == 0) {
							my $freq = sprintf("%.3f", $res_count{$k}/$clade1_count);
							$c1_residues = $c1_residues.$k."-".($freq*100)."% ";
						}
						else {
							$c1_residues = $c1_residues.$k."-".$res_count{$k}." ";
						}
					}
					$c1_residues =~ s/\s+$//;
					if ($only_count == 0) {
						print "Position $c2_cons_pos: $res2 (clade1: $c1_residues)\n"; 
						print LOG "Position $c2_cons_pos: $res2 (clade1: $c1_residues)\n";
					}
					$c2_spec_cons_count++;
					last;			
				}
				my $c1_cons_pos = $univ_cons_short1[$i2];
				if ($c2_cons_pos == $c1_cons_pos) {
					last;
				}
			}
		}
		print "Positions conserved in the second clade (type-I sites): $c2_spec_cons_count\n";
		print LOG "Positions conserved in the second clade (type-I sites): $c2_spec_cons_count\n";
	}
	
	print "\nAnalysis finished.\n";
	print LOG "\nAnalysis finished.\n";
	close LOG;
}
else {
	print "\nNo file specified or file does not exists. Pass a filename after the \"-f\" parameter.\n";
}

exit;

sub get_univ_cons {
	#TRANSFORM SEQUENCES INTO STRINGS
	my $seq1 = $_[0];
	my $seq2 = $_[1];
	my $arr = [];
	my $str1;
	my $str2;
	my $seq_str1 = IO::String->new(\$str1);
	my $seq_str2 = IO::String->new(\$str2);
	my $seqOut1 = Bio::SeqIO->new(-format => 'fasta', -fh => $seq_str1 );
	my $seqOut2 = Bio::SeqIO->new(-format => 'fasta', -fh => $seq_str2 );
	$seqOut1->write_seq($seq1);
	$seqOut2->write_seq($seq2);

	#REMOVE THE FASTA HEADER AND WHITESPACES, CONVERT TILDES TO HYPHENS
	$str1 = substr($str1,index($str1,"\n")+1);
	$str1 =~ s/[\n\r\s]+//g;
	$str1 =~ s/~/-/g;

	$str2 = substr($str2,index($str2,"\n")+1);
	$str2 =~ s/[\n\r\s]+//g;
	$str2 =~ s/~/-/g;	
			
	my @s1 = split(//,$str1);
	my @s2 = split(//,$str2);
	my $l = length ($str1);

	#COUNT IDENTITIES (EXCLUDE COMMON GAPS, AND COUNT UNKNOWNS)
	for (my $i = 0; $i < $l;$i++) {
		if ($s1[$i] eq $s2[$i] or ($s1[$i] eq "X") or ($s2[$i] eq "X")) {
			if ($s1[$i] ne "-") {
				push(@$arr,($i+1));
			}
		}
	}
	return $arr;
}

sub check_EI_pos {
	my @EI = ("ST","VI","SA","NS","DE","IL","NT","YF","EQ","LM","TA","RK","KQ","NH","GA","QP","SG","QH","VL","RH","AP","KN","RQ","SP","AV","DN","TM","TP","KT","VM","EA","SC","RS","RT","IM","QL","LW","PH","TI","LF","SL","KI","HY","DA","DH","LH","KM","RP","EG","VF","EK","DG","IF","SI","GV","RG","EV","SY","RI","RM","RL","GC","PL","RC","NY","SW","SF","DV","CF","NI","CW","CY","RW","GW","DY");
	my $sub = $_[0];
	my $ei_pos = 0;
	for (my $i1 = 0; $i1 < 75; $i1++) {
		if ($sub eq $EI[$i1] or reverse($sub) eq $EI[$i1]) {
			$ei_pos = $i1+1;
			last;
		}
	}
	return $ei_pos;
}



