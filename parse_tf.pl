#!/usr/bin/perl -w

use strict;

my $input = "";
my $out_path = "";
my $strand = "";
my $algo = "yeastract";
my $coreg = 0;
my $distrib = 0;
my $treshold = 0;
my $split = 0;
my $tresh1 = 0;
my $tresh2 = 0;
my $print_skip = 0;
my $distrib_short = 0;
my $distrib_table = 0;

if ($#ARGV == -1) {
	print "\n";
	print "This script parses transcription factor binding data from various sources and summarizes putative co-regulation statistics.\n";
	print "\nSyntax:\n";
	print " -f X		use file X as input\n";
	print " -out X		set a custom outfile\n";
	print " -rev		check the reverse strand\n";
	print " -coreg		check for coregulation in the data\n";
	print " -X 		set the data source: 'yeastract' for Yeastract tabular format (default),'ali' for AliBaba2 (Grabe, 2002) and 'promo' for Promo (Messeguer, 2002), 'explain' for Explain\n";
	print " -distrib	print TF distribution across all analyzed sequences\n";
	print " -treshold X_Y_Z	split the source data into 2 groups, starting at sequence no. X. Count only if at least Y and Z copies are present in the 1st and second group, respectively (use with -distrib)\n";
	print " -distrib_short	print just the names of identified TFs (use with -treshold)\n";
	print " -print_skip	print the skipped TFs (use with -treshold)\n";
	print " -distrib_table	print the TF distribution in a table\n";	
	print "\n";
	exit;
}

foreach my $argnum (0 .. $#ARGV) {
	if ($ARGV[$argnum] eq "-f") {
		$input = $ARGV[$argnum+1];
		if (-e $input) {
			print "\nLoading file $input...\n";
		}
		else {
			$input = "";
		}
	}
	if ($ARGV[$argnum] eq "-out") {
		$out_path = $ARGV[$argnum+1];
	}
	if ($ARGV[$argnum] eq "-rev") {
		$strand = "R";
	}
	if ($ARGV[$argnum] eq "-coreg") {
		$coreg = 1;
	}
	if ($ARGV[$argnum] eq "-ali") {
		$algo = "ali";
	}
	if ($ARGV[$argnum] eq "-promo") {
		$algo = "promo";
	}
	if ($ARGV[$argnum] eq "-patch") {
		$algo = "patch";
	}
	if ($ARGV[$argnum] eq "-jaspar") {
		$algo = "jaspar";
	}
	if ($ARGV[$argnum] eq "-explain") {
		$algo = "explain";
	}
	if ($ARGV[$argnum] eq "-distrib") {
		$distrib = 1;
	}
	if ($ARGV[$argnum] eq "-treshold") {
		$treshold = 1;
		my @tmp = split "_", $ARGV[$argnum+1];
		$split = $tmp[0];
		$tresh1 = $tmp[1];
		$tresh2 = $tmp[2];
	}
	if ($ARGV[$argnum] eq "-print_skip") {
		$print_skip = 1;
	}
	if ($ARGV[$argnum] eq "-distrib_short") {
		$distrib = 1;
		$distrib_short = 1;
	}
	if ($ARGV[$argnum] eq "-distrib_table") {
		$distrib = 1;
		$distrib_table = 1;
	}	
}

if ($input) {			
	my @seqs = ();
	my %TFs = ();
	my $target = "";
	my $input_short = substr($input, 0, -4);
	my $str = "";
	if ($distrib_short == 1) {
		$str = "_short";
	}
	$out_path = $input_short."_".$out_path.$str."_tf.txt";
	open (LOG, ">$out_path");


	#LOCATE SEQUENCE NAMES
	print "\nParsing input file $input for sequence names...\n";
	print LOG "\nParsing input file $input for sequence names...\n";

	open (TF_DATA, $input);	
	while (<TF_DATA>) {	
		my @fields = split(' ');	
		if ($algo eq "yeastract") {
			if (@fields == 5 and $fields[0] eq "Target") {
				push (@seqs, $fields[2]);
			}
		}
		elsif ($algo eq "promo") {
			if (@fields == 1 and substr($fields[0],0,1) eq ">") {
				push (@seqs, substr($fields[0],1));
			}
		}
		elsif ($algo eq "patch") {
			if (@fields > 2 and $fields[0] eq "Scanning") {
				push (@seqs, $fields[2]);
			}		
		}
		elsif ($algo eq "jaspar") {
			if (@fields > 2 and $fields[1] eq "putative") {
				push (@seqs, $fields[12]);
			}		
		}
		elsif ($algo eq "ali") {
			if (@fields == 2 and $fields[0] eq "Sequence") {
				push (@seqs, $fields[1]);
			}		
		}		
		elsif ($algo eq "explain") {
			if (@fields == 3 and $fields[0] eq "Gene") {
				push (@seqs, $fields[2]);
			}		
		}		
	}	
	close (TF_DATA);
		
	#LOCATE TF BINDING SITES
	print "\nParsing input file $input for TF binding data...\n";
	print LOG "\nParsing input file $input for TF binding data...\n";

	open (TF_DATA, $input);
	my $current_seq = "";
	while (<TF_DATA>) {
		my @fields = split(' ');		

		if ($algo eq "promo") {
			if (@fields == 9 and $fields[5] and $fields[1] ne "Factors") {
				my $seq = substr($fields[0],0,-1);
				push @{$TFs{$fields[1]}}, ("$seq\(".substr($fields[3],0,-1)."\)");		
			}
		}	
		elsif ($algo eq "yeastract") {
			if (@fields == 5 and $fields[0] eq "Target") {
				$current_seq = $fields[2];
			}
			if (@fields >= 4 and $current_seq ne "" and $fields[0] ne "Target" and $fields[0] ne "Back" and ($strand eq "" or $fields[3] eq $strand)) {	
				my $tf_count = @fields-4;
				push @{$TFs{$fields[0]}}, ("$current_seq\($fields[(2+$tf_count)]\)");
			}
		}	
		elsif ($algo eq "patch") {
			if (@fields == 3 and $fields[0] eq "Scanning") {
				$current_seq = $fields[2];
			}
			if ($current_seq ne "" and @fields == 7 and $fields[4] eq "100.00") {
				push @{$TFs{$fields[5]}}, ("$current_seq\($fields[1]\)");
			}
		}
		elsif ($algo eq "jaspar") {
			if (@fields == 13 and $fields[1] eq "putative") {
				$current_seq = $fields[12];
			}
			if ($current_seq ne "" and @fields == 8) {
				push @{$TFs{$fields[1]}}, ("$current_seq\($fields[4]\)");
			}
		}
		elsif ($algo eq "ali") {
			if (@fields == 2 and $fields[0] eq "Sequence") {
				$current_seq = $fields[1];
			}		
			if (@fields == 4 and $fields[0] ne "Class") {
				push @{$TFs{$fields[1]}}, ("$current_seq\($fields[2]\)");
			}		
		}
		elsif ($algo eq "explain") {
			if (@fields == 3 and $fields[0] eq "Gene") {
				$current_seq = $fields[2];
			}
			if (@fields == 6 and $fields[0] ne "(strand)" and $fields[2] >= 0.9) {
				push @{$TFs{$fields[5]}}, ("$current_seq\($fields[1]\)");
			}			
		}			
	}
	close (TF_DATA);

	#PRINT TF BINDING SITE DISTRIBUTION
	
	my ($TF, $sp);
	my $skipped = "";
	my $short = "";
	my $skip_count = 0;
	

	my @table = ();
	if ($distrib_table == 1) {
		$table[0] = "";
		for (my $i1 = 0; $i1 < @seqs; $i1++) {
			$table[$i1+1] = $seqs[$i1];
		}
		
	}
	while (($TF, $sp) = each(%TFs)) {
		if ($distrib == 1) {

			my $tfbs = "";
			my $pass = 0;
			foreach my $tf (@{$sp}) {
				$tfbs = $tfbs.$tf;
			}
			if ($treshold == 1) {
				my $count1 = 0;
				my $count2 = 0;
				for (my $i1 = 0; $i1 < @seqs; $i1++) {
					if ($tfbs =~ $seqs[$i1]) {
						if ($i1 < $split) {
							$count1++;
						}
						else {
							$count2++;
						}
					}
					if ($count1 >= $tresh1 and $count2 >= $tresh2) {	
						$pass = 1;
						$short = $short.$TF." ";
#							print $count1."__".$count2."\/";
						last;
					}
				}

			}	
			if ($distrib_short == 0 and ($treshold == 0 or $pass == 1)) {
				print "\n$TF\:\n";
				print LOG "\n$TF\:\n";
				if ($distrib_table == 1) {
					$table[0] = $table[0]."\t$TF";	
				}
				for (my $i1 = 0; $i1 < @seqs; $i1++) {				
					if ($tfbs =~ $seqs[$i1]) {
						if ($distrib_table == 1) {
							$table[$i1+1] = $table[$i1+1]."\t+";
						}
						print "$seqs[$i1] +\n";				
						print LOG "$seqs[$i1] +\n";				
					}
					else {
						if ($distrib_table == 1) {
							$table[$i1+1] = $table[$i1+1]."\t-";					
						}
						print "$seqs[$i1] -\n";				
						print LOG "$seqs[$i1] -\n";				
					}
				}
			}
			if ($treshold == 1 and $pass == 0) {
				$skipped = $skipped.$TF." ";
				$skip_count++;
			}
		}
		else {
			print "\n$TF:\t @{$sp}\n";
			print LOG "\n$TF:\t @{$sp}\n";
		}
	}
	if ($distrib_short == 1) {
		print "\nTFs identified by splitting at sequence number $split ($seqs[$split]), and using treshold $tresh1 before the split and $tresh2 after the split:\n\n$short\n\n";
		print LOG "\nTFs identified by splitting at sequence number $split ($seqs[$split]), and using treshold $tresh1 before the split and $tresh2 after the split:\n\n$short\n\n";
	}

	if ($print_skip == 1) {
		print "\n$skip_count TFs skipped by splitting at sequence number $split ($seqs[$split]), and using treshold $tresh1 before the split and $tresh2 after the split:\n\n$skipped\n\n";
		print LOG "\n$skip_count TFs skipped by splitting at sequence number $split ($seqs[$split]), and using treshold $tresh1 before the split and $tresh2 after the split:\n\n$skipped\n\n";

	}
	
	if ($distrib_table == 1) {
		print LOG "\n\n";		
		for (my $i1 = 0; $i1 < @table; $i1++) {
			print LOG $table[$i1]."\n";
		}
		print LOG "\n\n";
	}	


	my $seq_count = @seqs;
	if ($coreg == 1) {
		my @coreg_info;
		for (my $i1=0;$i1<$seq_count;$i1++) {
			for (my $i2=$i1+1;$i2<$seq_count and $i1!=$i2;$i2++) {
				my $coreg_count = 0;
				my ($TF2, $sp2);
				while (($TF2, $sp2) = each(%TFs)) {
					my $TF_tmp = join(" ",@{$sp2});
					if ($TF_tmp =~ $seqs[$i1] and $TF_tmp =~ $seqs[$i2]) {
						$coreg_count++;
					}
				}
				my $s = "$seqs[$i1] and $seqs[$i2] are coregulated by $coreg_count factors\n";
				push (@coreg_info, $s);
				print $s;
				print LOG $s;
			}
		}
		print "\n";
		print LOG "\n";
		foreach my $seq (@seqs) {
			my $c = 0;
			foreach my $s (@coreg_info) {
				my @f = split (' ', $s);
				if ($f[0] =~ $seq or $f[2] =~ $seq) {
					$c += $f[6];
				}
			}
			my $avg = $c/$seq_count;
			print "$seq is regulated on average by $avg factors\n";
			print LOG "$seq is regulated on average by $avg factors\n";
		}
	}
	print "\n";
	print LOG "\n";
	close LOG;
	close TF_DATA;
}
else {
	print "\nNo file specified or file does not exists. Pass a filename after the \"-f\" parameter.\n";
}


exit;


