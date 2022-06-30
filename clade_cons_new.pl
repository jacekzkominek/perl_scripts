#!/usr/bin/perl

use strict;
use warnings;
use IO::String;
use Bio::Seq; 
use Bio::SeqIO;
use Bio::AlignIO;
use Cwd 'abs_path';

my $input = "";
my $clade_split = 1;
my $align_split_pos = 0;
my $seq_split_pos = 0;
my $seq_split_ref = 1;
my $print_pos = 0;
my $print_id_pos = 0;
my $tolerance = 0;
my $matrix = 0;
my $diff_res_t1 = 0;
my $blosum = 62;
my $script_path = abs_path($0);
my $refseq = 0;
my $check_gaps = 0;

if ($#ARGV == -1) {
	print "\n";
	print "This script identifies conserved positions in an alignment. Positions conserved with clade-specific conservation between 2 groups of sequences.\n";
	print "\nSyntax:\n";
	print " -f X		use file X as input\n";
	print " -clade_split X	use all sequences before X-th as the first group, and the rest (including X-th) as the second group.\n";
#!!!	print " -align_split X  separately calculate conserved positions before and after position X in the alignment.\n";
#!!!	print " -seq_split X_Y	separately calculate conserved positions before and after position X in sequence Y (if Y is skipped, first sequence is taken.\n";
	print " -print_pos	print the individual positions of the main 3 classes (c1,c2 and c1c2_diff).\n";
	print " -print_id_pos	print the individual positions identical in both clades.\n";
	print " -tolerance X	allow X non-conforming residues in the annotation process.\n";
	print " -matrix 	create a position-specific matrix with class assignment.\n";
	print " -diff_res_t1 	require a different dominant residue for assignment to type-I site.\n";
	print " -blosum X	get substitution scores from a specific BLOSUM matrix (30, 35..[62]..100).\n";
	print " -refseq 	skip the first sequence and use its positions instead of alignment positions.\n";
	print " -check_gaps 	Do not skip gaps, treat them as 21st residue type.\n";
	print "\n";
	exit;
}

#LOAD PARAMETERS
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
	if ($ARGV[$argnum] eq "-print_pos") {
		$print_pos = 1;
	}
	if ($ARGV[$argnum] eq "-print_id_pos") {
		$print_id_pos = 1;
	}
	if ($ARGV[$argnum] eq "-tolerance") {
		$tolerance = $ARGV[$argnum+1];
	}
	if ($ARGV[$argnum] eq "-matrix") {
		$matrix = 1;
	}
	if ($ARGV[$argnum] eq "-diff_res_t1") {
		$diff_res_t1 = 1;
	}
	if ($ARGV[$argnum] eq "-blosum") {
		$blosum = $ARGV[$argnum+1];
	}
	if ($ARGV[$argnum] eq "-refseq") {
		$refseq = 1;
	}
	if ($ARGV[$argnum] eq "-check_gaps") {
		$check_gaps = 1;
	}
}

if ($input and $clade_split > 1) {
	open (LOG, ">$input\_log");

	print "\nLoading sequences from $input...\n";
	print LOG "\nLoading sequences from $input...\n";
	
	#READ IN THE BLOSUM MATRIX
	my @blosum_matrix = ();
	print "\nLoading the blosum$blosum matrix...\n";
	my $script_dir = substr($script_path, 0, rindex($script_path,"\/"));
	print "";
	open BLOSUM, "$script_dir\/blosum/blosum$blosum\.blo" or die $!;
	while (<BLOSUM>) {
		if (!($_ =~ "#")) {
			push @blosum_matrix, $_;	
		}
	}
	
	my $seq_in = Bio::AlignIO->new(-format => 'fasta', -file => $input);
	my $seq_aln = $seq_in->next_aln;
	my $seq_count = $seq_aln->num_sequences;
	my $aln_len = $seq_aln->length;
	if ($seq_count < 2) {
		print "\nLess than 2 sequences in the file. Exiting\n";
		exit;
	}
	my @matrix = ();
	
	#SPLIT THE SEQUENCES INTO 2 CLADES
	my $clade1 = new Bio::SimpleAlign();
	my $clade2 = new Bio::SimpleAlign();
	my $start_seq = 1;
	if ($refseq == 1) {
		$start_seq = 2;
	}
	for (my $i1 = $start_seq; $i1 <= $seq_count; $i1++) {
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
	
	#PRINT SEQUENCE TITLES
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

	my @pos_c1_cons = ();
	my @pos_c1_cons_res = ();
	my @pos_c2_cons = ();
	my @pos_c2_cons_res = ();
	my @pos_c1c2_cons = ();
	my @pos_c1c2_cons_res = ();
	my @pos_c1c2_diff = ();
	my @pos_c1c2_diff_res = ();

	#PARSE ALL POSITIONS AND ASSIGN THEM TO CLASSES
	for (my $pos = 1; $pos <= $aln_len; $pos++) {
		my %c1_res_count = ();
		my %c2_res_count = ();
		my @sorted_c1_res_count = ();
		my @sorted_c1_res_type = ();
		my @sorted_c2_res_count = ();
		my @sorted_c2_res_type = ();

		
		#PARSE CLADE1 AND STORE RESIDUE COUNTS. SORT.
		for (my $i1 = 1; $i1 <= $clade1_count; $i1++) {
			my $seq1 = $clade1->get_seq_by_pos($i1);
			my $seq1_res = $seq1->subseq($pos,$pos);
			if (!(defined $c1_res_count{$seq1_res})) {
				$c1_res_count{$seq1_res} = 1;
			}
			elsif ($c1_res_count{$seq1_res} > 0) {
				$c1_res_count{$seq1_res} = $c1_res_count{$seq1_res}+1;
			}	
		}
		
		if ($check_gaps == 0 and defined $c1_res_count{"-"}) {
			next;
		}

		foreach my $key (sort {$c1_res_count{$b} <=> $c1_res_count{$a}} keys %c1_res_count) {
			push @sorted_c1_res_count, $c1_res_count{$key};
			push @sorted_c1_res_type, $key;
		}
		#PARSE CLADE2 AND STORE RESIDUE COUNTS. SORT.
		for (my $i1 = 1; $i1 <= $clade2_count; $i1++) {
			my $seq2 = $clade2->get_seq_by_pos($i1);
			my $seq2_res = $seq2->subseq($pos,$pos);
			if (!(defined $c2_res_count{$seq2_res})) {
				$c2_res_count{$seq2_res} = 1;
			}
			elsif ($c2_res_count{$seq2_res} > 0) {
				$c2_res_count{$seq2_res} = $c2_res_count{$seq2_res}+1;
			}
		}
		
		if ($check_gaps == 0 and defined $c2_res_count{"-"}) {
			next;
		}		
		foreach my $key (sort {$c2_res_count{$b} <=> $c2_res_count{$a}} keys %c2_res_count) {
			push @sorted_c2_res_count, $c2_res_count{$key};
			push @sorted_c2_res_type, $key;
		}			

		my $matrix_line = "$pos";
		#ASSIGN POSITIONS TO CLASSES BASED ON CONSERVATION (NO CONSERVATION = NO ASSIGNMENT)
		if ($sorted_c1_res_count[0] >= $clade1_count-$tolerance and $sorted_c2_res_count[0] >= $clade2_count-$tolerance) {
			#CLASS 1: CONSERVED AND IDENTICAL IN BOTH CLADES
			if ($sorted_c1_res_type[0] eq $sorted_c2_res_type[0]) {
				push @pos_c1c2_cons, $pos;
									
				my $c1c2_cons_res = "\tclade1: ";
				for (my $i1=0; $i1 < @sorted_c1_res_count; $i1++) {
					$c1c2_cons_res = $c1c2_cons_res.$sorted_c1_res_type[$i1]."-".$sorted_c1_res_count[$i1]." ";
				}
				$c1c2_cons_res = $c1c2_cons_res."\tclade2: ";					
				for (my $i1=0; $i1 < @sorted_c2_res_count; $i1++) {
					$c1c2_cons_res = $c1c2_cons_res.$sorted_c2_res_type[$i1]."-".$sorted_c2_res_count[$i1]." ";
				}
				push @pos_c1c2_cons_res, $c1c2_cons_res;
				$matrix_line = $matrix_line."\t1\t0\t0\t0\t0";
			}
			#CLASS 2: CONSERVED BUT DIFFERENT IN BOTH CLADES				
			else {
				push @pos_c1c2_diff, $pos;
				
				
				my $c1c2_diff_res = "clade1: ";
				for (my $i1=0; $i1 < @sorted_c1_res_count; $i1++) {
					$c1c2_diff_res = $c1c2_diff_res.$sorted_c1_res_type[$i1]."-".$sorted_c1_res_count[$i1]." ";
				}
				$c1c2_diff_res = rtrim($c1c2_diff_res);
				$c1c2_diff_res = $c1c2_diff_res."\tclade2: ";					
				for (my $i1=0; $i1 < @sorted_c2_res_count; $i1++) {
					$c1c2_diff_res = $c1c2_diff_res.$sorted_c2_res_type[$i1]."-".$sorted_c2_res_count[$i1]." ";
				}
				$c1c2_diff_res = rtrim($c1c2_diff_res);				
				#CHECK EVOLUTIONARY INDEX
				my $ei_rank = "";
				if (check_EI_rank($sorted_c1_res_type[0].$sorted_c2_res_type[0]) > 0) {
					$ei_rank = check_EI_rank($sorted_c1_res_type[0].$sorted_c2_res_type[0]);
				}
				else {
					$ei_rank = "XX";
				}

				#CHECK BLOSUM SCORE
				my $blosum_score = check_blosum_score(\@blosum_matrix,$sorted_c1_res_type[0],$sorted_c2_res_type[0]);
				
				push @pos_c1c2_diff_res, "\t".$c1c2_diff_res."\t\t(".$sorted_c1_res_type[0]."<->".$sorted_c2_res_type[0]."\tEI:$ei_rank\tB$blosum:$blosum_score)";
				$matrix_line = $matrix_line."\t0\t0\t0\t1\t0";
			}
		}		
		#CLASS 3: CONSERVED ONLY IN THE FIRST CLADE
		elsif ($sorted_c1_res_count[0] >= $clade1_count-$tolerance) {
			if ($diff_res_t1 == 0 or ($diff_res_t1 == 1 and $sorted_c1_res_type[0] ne $sorted_c2_res_type[0])) {
				push @pos_c1_cons, $pos;
				#GET BLOSUM AND EI SCORES FOR MOST FREQUENT RESIDUES
				my $blosum_score = check_blosum_score(\@blosum_matrix,$sorted_c1_res_type[0],$sorted_c2_res_type[0]);
				my $ei_rank = "";
				if (check_EI_rank($sorted_c1_res_type[0].$sorted_c2_res_type[0]) > 0) {
					$ei_rank = check_EI_rank($sorted_c1_res_type[0].$sorted_c2_res_type[0]);
				}
				else {
					$ei_rank = "XX";
				}	
				
				my $c1_cons_res = "\tclade1: ";
				for (my $i1=0; $i1 < @sorted_c1_res_count; $i1++) {
					$c1_cons_res = $c1_cons_res.$sorted_c1_res_type[$i1]."-".$sorted_c1_res_count[$i1]." ";
				}
				$c1_cons_res = rtrim($c1_cons_res);
				$c1_cons_res = $c1_cons_res."\tclade2: ";
				for (my $i1=0; $i1 < @sorted_c2_res_count; $i1++) {
					$c1_cons_res = $c1_cons_res.$sorted_c2_res_type[$i1]."-".$sorted_c2_res_count[$i1]." ";
				}
				$c1_cons_res = rtrim($c1_cons_res);				
				$c1_cons_res = $c1_cons_res."\t\t(".$sorted_c1_res_type[0]."<->".$sorted_c2_res_type[0]."\tEI:$ei_rank\tB$blosum:$blosum_score)";
				push @pos_c1_cons_res, $c1_cons_res;
				$matrix_line = $matrix_line."\t0\t1\t0\t0\t0";	
			}		
		}
		#CLASS 4: CONSERVED ONLY IN THE SECOND CLADE			
		elsif ($sorted_c2_res_count[0] >= $clade2_count-$tolerance) {
			if ($diff_res_t1 == 0 or ($diff_res_t1 == 1 and $sorted_c2_res_type[0] ne $sorted_c1_res_type[0])) {
				push @pos_c2_cons, $pos;
				#GET BLOSUM AND EI SCORES FOR MOST FREQUENT RESIDUES
				my $blosum_score = check_blosum_score(\@blosum_matrix,$sorted_c1_res_type[0],$sorted_c2_res_type[0]);
				my $ei_rank = "";
				if (check_EI_rank($sorted_c1_res_type[0].$sorted_c2_res_type[0]) > 0) {
					$ei_rank = check_EI_rank($sorted_c1_res_type[0].$sorted_c2_res_type[0]);
				}
				else {
					$ei_rank = "XX";
				}	
				
				my $c2_cons_res = "\tclade2: ";
				for (my $i1=0; $i1 < @sorted_c2_res_count; $i1++) {
					$c2_cons_res = $c2_cons_res.$sorted_c2_res_type[$i1]."-".$sorted_c2_res_count[$i1]." ";
				}
				$c2_cons_res = rtrim($c2_cons_res);
				$c2_cons_res = $c2_cons_res."\tclade1: ";
				for (my $i1=0; $i1 < @sorted_c1_res_count; $i1++) {
					$c2_cons_res = $c2_cons_res.$sorted_c1_res_type[$i1]."-".$sorted_c1_res_count[$i1]." ";
				}
				$c2_cons_res = rtrim($c2_cons_res);				
				$c2_cons_res = $c2_cons_res."\t\t(".$sorted_c2_res_type[0]."<->".$sorted_c1_res_type[0]."\tEI:$ei_rank\tB$blosum:$blosum_score)";				
				push @pos_c2_cons_res, $c2_cons_res;
				$matrix_line = $matrix_line."\t0\t0\t1\t0\t0";	
			}
		}
		else {
			$matrix_line = $matrix_line."\t0\t0\t0\t0\t1";
		}
		push @matrix, $matrix_line;
	}
	my $c1c2_cons_size = @pos_c1c2_cons;	
	my $c1c2_diff_size = @pos_c1c2_diff;
	my $c1_cons_size = @pos_c1_cons;
	my $c2_cons_size = @pos_c2_cons;

	print "\nPositions conserved in both clades: $c1c2_cons_size \/ $aln_len \n";
	print LOG "\nPositions conserved in both clades: $c1c2_cons_size \/ $aln_len \n";
	if ($print_id_pos == 1) {
		for (my $i1 = 0; $i1 < $c1c2_cons_size;$i1++) {
			my $pos = $pos_c1c2_cons[$i1];
			print "#".($i1+1)." Position\t$pos:$pos_c1c2_cons_res[$i1]\n";
			print LOG "#".($i1+1)." Position\t$pos:$pos_c1c2_cons_res[$i1]\n";
		}
	}

	print "\nPositions conserved but different in both clades: $c1c2_diff_size \/ $aln_len \n";
	print LOG "\nPositions conserved but different in both clades: $c1c2_diff_size \/ $aln_len \n";		
	if ($print_pos == 1) {
		for (my $i1 = 0; $i1 < $c1c2_diff_size;$i1++) {
			my $pos = $pos_c1c2_diff[$i1];			
			print "#".($i1+1)." Position\t$pos:$pos_c1c2_diff_res[$i1]\n";
			print LOG "#".($i1+1)." Position\t$pos:$pos_c1c2_diff_res[$i1]\n";
		}
	}
	print "\nPositions conserved only in the first clade: $c1_cons_size \/ $aln_len \n";
	print LOG "\nPositions conserved only in the first clade: $c1_cons_size \/ $aln_len \n";
	if ($print_pos == 1) {
		for (my $i1 = 0; $i1 < $c1_cons_size;$i1++) {
			my $pos = $pos_c1_cons[$i1];			
			print "#".($i1+1)." Position\t$pos:$pos_c1_cons_res[$i1]\n";
			print LOG "#".($i1+1). " Position\t$pos:$pos_c1_cons_res[$i1]\n";
		}
	}
	print "\nPositions conserved only in the second clade: $c2_cons_size \/ $aln_len \n";
	print LOG "\nPositions conserved only in the second clade: $c2_cons_size \/ $aln_len \n";
	if ($print_pos == 1) {
		for (my $i1 = 0; $i1 < $c2_cons_size;$i1++) {
			my $pos = $pos_c2_cons[$i1];			
			print "#".($i1+1)." Position\t$pos:$pos_c2_cons_res[$i1]\n";
			print LOG "#".($i1+1)." Position\t$pos:$pos_c2_cons_res[$i1]\n";
		}
	}
	if ($matrix == 1) {
		open (MATRIX, ">$input\_class_matrix.txt");
		print MATRIX "Pos\tT0\tT1A\tT1B\tT2\tN\n";
		foreach my $line (@matrix) {
			print MATRIX $line."\n";
		}
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

sub check_EI_rank {
	my @EI = ("ST","VI","SA","NS","DE","IL","NT","YF","EQ","LM","TA","RK","KQ","NH","GA","QP","SG","QH","VL","RH","AP","KN","RQ","SP","AV","DN","TM","TP","KT","VM","EA","SC","RS","RT","IM","QL","LW","PH","TI","LF","SL","KI","HY","DA","DH","LH","KM","RP","EG","VF","EK","DG","IF","SI","GV","RG","EV","SY","RI","RM","RL","GC","PL","RC","NY","SW","SF","DV","CF","NI","CW","CY","RW","GW","DY");
	my $sub = $_[0];
	my $ei_rank = 0;
	for (my $i1 = 0; $i1 < 75; $i1++) {
		if ($sub eq $EI[$i1] or reverse($sub) eq $EI[$i1]) {
			$ei_rank = $i1+1;
			last;
		}
	}
	return $ei_rank;
}


sub check_blosum_score {
	my $blosum_ref = $_[0];
	my @blosum_matrix = @$blosum_ref;
	my $res1 = $_[1];
	my $res2 = $_[2];
	
	my $blosum_score = "";
	my @aa_list = split(' ',$blosum_matrix[0]);
	my $col = "";
	for (my $i1 = 0; $i1 < @aa_list; $i1++) {
		if ($aa_list[$i1] eq $res1) {
			$col = $i1;
			last;
		}
	}
	foreach my $aa (@blosum_matrix) {
		if (substr($aa, 0, 1) eq $res2) {
			my @aa_scores = split(' ',$aa);
			$blosum_score = $aa_scores[int($col)+1];
			last;
		}
	}
	return $blosum_score;
}

sub rtrim($)
{
	my $string = shift;
	$string =~ s/\s+$//;
	return $string;
}
