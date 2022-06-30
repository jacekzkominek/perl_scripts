#!/usr/bin/perl -w

use strict;
use IO::String;
use Bio::AlignIO;
use Bio::SeqIO;
use Bio::Tools::Run::StandAloneBlastPlus;
use Bio::SearchIO;

my $input = "";
my $db = "";
my $out_dir = ".";
my $no_margin = 0;
my $margin = 300;
my $best_only = 0;
my $no_extract = 0;

if ($#ARGV == -1) {
	print "\n";
	print "This script runs BLAST and extracts hits from a specific database.\n";
	print "\nSyntax:\n";
	print " -f X		use file X as input\n";
	print " -db X		seach in blast database X\n";
	print " -out_dir X	save output into directory X\n";
	print " -best_only 	extract only 1 best hit\n";
	print " -no_margin  	extract only the hit itself\n";
	print " -margin X 	extract the hit itself and X positions upstream and downstream (150 by default)\n";	
	print " -no_extract 	only run BLAST, do not extract hits.\n";	
	print "\n";
	exit;
}

foreach my $argnum (0 .. $#ARGV) {
	if ($ARGV[$argnum] =~ "-f") {
		$input = $ARGV[$argnum+1];
	}
	if ($ARGV[$argnum] =~ "-db") {
		$db = $ARGV[$argnum+1];
	}
	if ($ARGV[$argnum] =~ "-out_dir") {
		$out_dir = $ARGV[$argnum+1];
	}
	if ($ARGV[$argnum] =~ "-best_only") {
		$best_only = 1;
	}
	if ($ARGV[$argnum] =~ "-no_margin") {
		$no_margin = 1;
	}
	if ($ARGV[$argnum] =~ "-margin") {
		$margin = $ARGV[$argnum+1];
	}
	if ($ARGV[$argnum] =~ "-no_extract") {
		$no_extract = 1
	}
}

my $seq_source = $input;
my $seq_source_short = substr($input,0,length($input)-4);
my $seq_source_dir = substr($input,0,rindex($input,"\/")+1);

#LOAD BLAST QUERY (SOURCE) AND THE OUTPUT FILE
print "\n";
my $seq_input = Bio::SeqIO->new(-file=>"$seq_source", -format => 'Fasta');
my $seq_out;
if ($best_only == 1) {
	$seq_out = Bio::SeqIO->new(-file=>">$db\_best_blast\_seqs.fas", -format => 'Fasta');
}
elsif ($no_extract == 0) {
	$seq_out = Bio::SeqIO->new(-file=>">$db\_blast\_seqs.fas", -format => 'Fasta');
}

while (my $seq = $seq_input->next_seq) {
	#MAIN LOOP PART1 - SET THE APPROPRIATE BLAST ALGORITHM AND EXECUTE IT
	my $prog;
	if ($seq_source =~ "_nn") {
		$prog = "blastn";
		if ($db =~ "_aa") {
			$prog = "blastx";
		}
	}
	else {
		$prog = "tblastn";
		if ($db =~ "_aa") {
			$prog = "blastp";
		}		
	
	}

	my $seq_name = $seq->display_id;
	print "\nBlasting $db for $seq_name...";
	my $out2 = $out_dir;
	my $out3 = "$out2/blast_files";
	mkdir($out3);
#	my @params = (program => $prog, -database => $db, -outfile => "$out3/blast\_$seq_name\_$db.txt");
	my @params = (-db_name => $db);
	my $outfile = "$out3/blast\_$seq_name\_$db.txt";
#	my $factory = Bio::Tools::Run::StandAloneBlast->new(@params);
	my $factory = Bio::Tools::Run::StandAloneBlastPlus->new(@params);
#	$factory->F("false");
#	$factory->e("1.0");
#	my $rep = $factory->blastall($seq);
	my $rep;
	if ($prog eq "blastp") {
		$rep = $factory->blastp(-query => $seq_source, -outfile => $outfile, -method_args => [-evalue => 1E-10]);
	}
	if ($prog eq "blastn") {
		$rep = $factory->blastn(-query => $seq_source, -outfile => $outfile, -method_args => [-evalue => 1E-10]);
	}	
	if ($prog eq "blastx") {
		$rep = $factory->blastx(-query => $seq_source, -outfile => $outfile, -method_args => [-evalue => 1E-10]);
	}	
	if ($prog eq "tblastn") {
		$rep = $factory->tblastn(-query => $seq_source, -outfile => $outfile, -method_args => [-evalue => 1E-10]);
	}			
	#MAIN LOOP PART2 - RETRIEVE THE (BEST) BLAST HITS SEQUENCES
	#BY PARSING THROUGH RESULT->HIT->HSP
	
	if ($no_extract == 0) {
		if ($best_only == 1) {
			print "finished\nRetrieving only best hit sequences...";
		}
		else {
			print "finished\nRetrieving hit sequences...";
		}
#		my $res = $rep->next_result();
		my $res = $factory->next_result();
		my $first_hit = 0;
		my $hit_num = 1;
		while (my $hit = $res->next_hit()) {
			my $n = $hit->name();
	#		my $hit_s = $hit->start('Sbjct');
	#		my $hit_e = $hit->end('Sbjct');

			#ONLY PROCEED IF THE HIT IS OF LOGICAL SIZE (5000 FOR HSP70)
			#OR IT IS THE VERY FIRST HIT
	#		if ($first_hit == 0 or $hit_e-$hit_s < 5000) {	
				my $hsp_num = 1;
				while (my $hsp = $hit->next_hsp) {
					my $hsp_s = $hsp->start('Sbjct');
					my $hsp_e = $hsp->end('Sbjct');
					my $hit_bit_score = $hsp->bits();
					my $len = $hsp->length('Sbjct');

					#IF THE SEQUENCE WILL BE RETRIEVED WITHOUT
					#ANY SEQUENCE MARGIN THEN TAKE IT FROM THE HIT ITSELF
	#				if ($no_margin == 1) {
	#					my $id = $hit->name();
	#					my $seq2 = $hsp->seq('hit');

	#					my $new_id =("$seq_name\_hit$hit_num\.$hsp_num\|Length\:$len\|Score\:$hit_bit_score\|$db\_$id\[$hsp_s-$hsp_e\]");
	#					$seq2->id($new_id);
	#					$seq_out->write_seq($seq2);
	#				}
	#				else{			
						#RETRIEVE THE SEQUENCE FROM THE DATABASE WITH OR WITHOUT A MARGIN
						my $db_file = "/media/sda4/work/db/blast/$db\.fas";
						my $seq_db = Bio::SeqIO->new(-file=>$db_file, -format => 'Fasta');
						while (my $seq2 = $seq_db->next_seq) {
							my $id = $seq2->id;
							if ($id eq $n) {
								my $fin = 0;
								if ($no_margin == 0) {
									if ($hsp_s <= $margin) {
										$hsp_s = 1;
									}
									else {
										$hsp_s = $hsp_s-$margin;
									}
									my $l = $seq2->length();
									if ($hsp_e+$margin > $l) {
										$hsp_e = $l;
										$fin = 1;
									}
									else {
										$hsp_e = $hsp_e+$margin;
									}
								}
								my $seq3 = $seq2->subseq($hsp_s,$hsp_e);
								if ($fin == 1) {
									$hsp_e = $hsp_e."E";
								}
								my $seq4 = Bio::Seq->new(-format => "Fasta",-seq => $seq3);
								my $new_id =("$seq_name\_hit\_$hit_num\.$hsp_num\|Length\:$len\|Score\:$hit_bit_score\|$db\_$id\[$hsp_s-$hsp_e\]");
								$seq4->id($new_id);
								$seq_out->write_seq($seq4);
							}
						}
	#				}
					$hsp_num++;
					if ($best_only == 1) {
						last;
					}
				}
				#SET A MARK THAT THE FIRST HIT IS PAST
				$first_hit = 1;
				$hit_num++;
	#		}
			if ($best_only == 1) {
				last;
			}
		}
		print "finished\n";
	}
}

exit;

