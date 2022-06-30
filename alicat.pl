#!/usr/bin/env perl

=head1 NAME

alicat.pl

=head1 SYNOPSIS

alicat.pl [-a] [-A] [-c chr] [-d] [-g] [-l loOffset] [-h hiOffset] [-i pats] [-o file] [-p] [-q qTh] [-r region] [-s] [-t len [-L] [-P]] [-R] [-v] [-x pats] [-y] [-z] file1 file2 ...

=head1 DESCRIPTION

Alignment file printing script, originally inspired by "cat -n". The
input file(s) are as produced by the SGRP pipeline, e.g. layout.gz,
sequenced.gz, imputed.gz. See ftp://ftp.sanger.ac.uk/pub/dmc/yeast/latest.

* Specifying the positions to be printed.

If -l and/or -h are specified, only positions within that region are
output.

See also -c below for specifying the chromosome or contig file.

If -a is given, offsets refer to positions in the alignment,
i.e. including reference-sequence gaps; otherwise (the default) gaps
do not add to the offset and we show such lines with offsets NNN.001,
NNN.002 etc. To handle insertions of length greater than 999 using a
three-character extension, we use A-Z for 1000 to 3599 and then a-z
for 3600 to 6199, and * for any higher values, so 1234 is C34 and 3634
is a34. (This is messy, but better than changing the layout, which
existing code may rely on).

If -A is given, positions are always printed, even when they consist
entirely of gaps. Otherwise, gap-only positions are not printed (but
positions that have only gaps and "N" symbols _are_ printed).  If you
want to be able to compare line-by-line alicat.pl output from two
parallel files, such as sequenced.gz and imputed.gz for the same
chromosome, you should use -A.

If -z is given, count first position as 0, else 1.

If -r is given, the value should be of the form
identifier:fileName.gff. fileName.gff will be searched for an object
with the specified identifier, the low and high ends of the first such
object found will be added to -l and -h respectively, and if the
object is on the reverse strand, -R (see below) will be set. Thus if
-l and -h are omitted, you get exactly the region covered by the
object in question, and if they are present, they adjust that region
(see examples below). The chromosome for the region will also be used
as a filter on the fileN (data) files in the command: only the (first)
file with a directory name component equal to the chromosome will be
used.

If -L is given, the input file(s) file1 file2 ... are treated as being
in "layout" format (see the SGRP manual). Otherwise, they are treated
this way only if their names contain the string "layout".

* Specifying the strains to be printed.

(These only work when there is one column pair per strain, e.g. "sorted"
and "imputed" files, not "layout" or "rqpw").

If -i str1,str2,... is given, include only strains matching one of those names
    (Perl-style regular expressions allowed) except any matching the -x list.
    The match must be to the WHOLE of the strain name, so if you want all
    the UWOPS strains, say -i 'UWOPS.*'

If -x str1,str2,... is given, exclude strains matching any of those names
(even if they are given in -i).

* Output destination

If -o is given, output to that file; default is standard output. If file name
ends in .gz, gzip to it.

* Position selection and gross line format

If -p is given, only polymorphic positions (those with more than one
different nucleotide of sufficient quality -- see -q) are printed.

If -q is given, replace all nucleotides with less than that quality
value by "N" (or "n" when lower case). If -p is also given, only count
values with at least that quality (as a number) when deciding
polymorphism.

If -s is given, separate nucleotides and qualities in output (see below).

If -v (but not -t, for which see below) is given, print strain names
vertically at the start, in the columns that will be used for their
sequences.

If -t len is given, turn the whole alignment to be horizontal, with
len alignment positions per line. Because each line will be prefixed by
a strain name, the line length will be greater than len.

If -P is given as well as -t, only print lines for strains that differ
from the reference.

* Modifying individual nucleotide and quality value characters

If -d is given, use ":" for a sequenced value that is the same as the
reference, and "." for an imputed value that is the same.

If -g is given, leave numeric gap characters (created during the 
imputation process) as they are; default is to translate them to '-'.

If -y is given, replace quality chars with the tens part of the quality value;
thus e.g. "P" is quality ord('P')-ord('!') = 47 so use 4. Leave space as space.

-c chr: If there is more than one input file, use only the first one
with a name component equal to chr. In addition, in vertical output mode,
each input line starting with ">" is printed as is, except that for
the first line, if -v is given, we print the strain names on it vertically.
Otherwise, the value of "-c" (if any) is put in after the
">". Thus if we say "alicat.pl -c foo bar", and the first ">" line is
">xyz", the output will be ">foo xyz". "-c" is set automatically
when "-r" succeeds (see above).

If -R is given, or is set implicitly by using "-r" with a reverse
strand object, nucleotide values are reverse complemented, and the
output lines are printed in reverse order. (Known problem, currently
not a priority to fix: this applies to the header line, giving the
strain name, too, so it occurs at the end of the file. Because of
this, -R doesn't work well with vertical output (-v or -t)).

* Vertical output format

The columns in vertical output (created without "-t") are as follows:

1-7:   the offset in the chromosome
8-11:  if the first character in the input line is a gap symbol
       ("-" or a digit), cols 8-11 are of the form ".NNN", and
        cols 1-7 are the same as in the previous line; otherwise they
        are blank. Thus MMM.NNN means "the NNN'th inserted position
        after offset MMM".
13-15: Consensus nucleotide value. If consensus is unanimously X,
       field is " X "; if not unanimous, i.e. if the line is polymorphic, 
       it's "[X]". A lower-case consensus means there is a gap in some 
       positions and nucleotides in others.
16:    always blank
17:    always "|"

For the remaining columns: N is the number of strains including
the reference, i.e. each non-">" input line has 2*N characters:

If "-s" ("separate") switch is given, then remaining columns are:

18 to 17+N: the nucleotide values from the input (with digits translated
            to "-" if "-t" was given).

18+N to 20+N: always "| |"
21+N to 20+2*N: the quality characters from the input.
21+2*N:       always "|"

If "-s" is not given, then:

18 to 17+2*N: an exact copy of the input, i.e. with nucleotides and qualities
              still interleaved (but digits are translated to "-" if "-t" was 
              given).
18+2*N:       always "|"

* Examples

1) alicat.pl -dy -v -t 100 -P -r GAL1:/foo/ref/cere/genome.gff p?/*/imputed_m.gz

Look for GAL1 in the specified GFF file, then find the first of the
imputed_m.gz files whose name contains (in the second field) the right
chromosome name. Convert and print that range of lines from that file.
Lines are turned (-t 100) to be horizontal, with 100 alignment
positions per line. Within each batch of turned lines, only print details
of strains that are different from the reference (-P) somewhere in the line.
Replace (-d) sequenced values identical to the reference by ":", and imputed
by ".". Replace (-y) each Phred quality character by its tens digit.

If /foo/ref/cere/genome.gff contains the line

chr02 SGD gene 279021 280607 . + . ID=YBR020W;Name=YBR020W;gene=GAL1;Alias=GAL1 ...

then the above command is equivalent each of

alicat.pl -dy -v -t 100 -P -l 279021 -h 280607 -c chr02 p?/*/imputed_m.gz
alicat.pl -dy -v -t 100 -P -l 279021 -h 280607 p2/chr02/imputed_m.gz

2) alicat.pl -r GAL1:/foo/ref/cere/genome.gff -l -2000 -h 200 p?/*/imputed_m.gz

Look for GAL1, and print information starting 2000 bp upstream from it and
ending 200 bp downstream.

=head1 AUTHOR

David Carter, C<dmc@sanger.ac.uk>. Suggestions for new functionality welcome.

=head1 TERMS AND CONDITIONS

Copyright (c) 2007 Wellcome Trust Sanger Institute (WTSI)

This program is free software; you can redistribute it and/or modify
it under the same terms as Perl itself. See the Artistic License file
in the main Perl distribution for specific terms and conditions of
use.  In addition, the following disclaimers apply:

WTSI makes no representations whatsoever as to the SOFTWARE contained
herein.  It is experimental in nature and is provided WITHOUT WARRANTY
OF MERCHANTABILITY OR FITNESS FOR A PARTICULAR PURPOSE OR ANY OTHER
WARRANTY, EXPRESS OR IMPLIED.  CSHL MAKES NO REPRESENTATION OR
WARRANTY THAT THE USE OF THIS SOFTWARE WILL NOT INFRINGE ANY PATENT OR
OTHER PROPRIETARY RIGHT.

By downloading this SOFTWARE, your Institution hereby indemnifies WTSI
against any loss, claim, damage or liability, of whatsoever kind or
nature, which may arise from your Institution's respective use,
handling or storage of the SOFTWARE.

If publications result from research using this SOFTWARE, we ask that
WTSI be acknowledged and/or credit be given to WTSI scientists, as
scientifically appropriate.

=head1 SEE ALSO

       

=cut

use strict;

use File::Basename;

use FileHandle;
use Getopt::Std;
use strict;

my %OPT;
getopts('aAc:dgG:i:Ll:h:o:pPq:r:Rst:vx:yzT:',\%OPT);

# Bit positions are:
#  A 1
#  C 2
#  G 4
#  T 8
#  - 16

my %TRANSLATION = ( TTT => 'Phe', TTC => 'Phe', TTA => 'Leu', TTG => 'Leu',
		    TCT => 'Ser', TCC => 'Ser', TCA => 'Ser', TCG => 'Ser', 
		    TAT => 'Tyr', TAC => 'Tyr', TAA => 'STP', TAG => 'STP', 
		    TGT => 'Cys', TGC => 'Cys', TGA => 'STP', TGG => 'Trp', 
		    CTT => 'Leu', CTC => 'Leu', CTA => 'Leu', CTG => 'Leu', 
		    CCT => 'Pro', CCC => 'Pro', CCA => 'Pro', CCG => 'Pro', 
		    CAT => 'His', CAC => 'His', CAA => 'Gln', CAG => 'Gln', 
		    CGT => 'Arg', CGC => 'Arg', CGA => 'Arg', CGG => 'Arg', 
		    ATT => 'Ile', ATC => 'Ile', ATA => 'Ile', ATG => 'Met',
		    ACT => 'Thr', ACC => 'Thr', ACA => 'Thr', ACG => 'Thr',
		    AAT => 'Asn', AAC => 'Asn', AAA => 'Lys', AAG => 'Lys', 
		    AGT => 'Ser', AGC => 'Ser', AGA => 'Arg', AGG => 'Arg',
		    GTT => 'Val', GTC => 'Val', GTA => 'Val', GTG => 'Val', 
		    GCT => 'Ala', GCC => 'Ala', GCA => 'Ala', GCG => 'Ala', 
		    GAT => 'Asp', GAC => 'Asp', GAA => 'Glu', GAG => 'Glu', 
		    GGT => 'Gly', GGC => 'Gly', GGA => 'Gly', GGG => 'Gly');

my @DNA_SYM;
my %DNA_BITS;

if (defined $OPT{T}) {
    # use as "aliturn" in pipe out of alicat. -T is deliberately not
    # documented, as you should not directly call it on the command line.
    aliturn($OPT{T},$OPT{P});
} else {
    alicat(\%OPT,@ARGV);
}

sub alicat {
    my ($pOpt,@files) = @_;
    setDNABitsAndSymbols();
    # Set the output handle: see openToWrite
    if ($pOpt->{r}) {
	# If we specified a region, go via alicatRegion to choose
	# the file to look in.
	alicatRegion($pOpt,@files);
    } else {
	my $gh = openToWrite($pOpt->{o},$pOpt->{t},$pOpt->{P},$pOpt->{R});
	foreach my $f (@files) {
	    # If -c CHR was specified, the file has to include CHR
	    # in its name.
	    if (!$pOpt->{c} || hasNameComponent($f,$pOpt->{c})) {
		# Do the work ...
		alicat1($f,$gh,$pOpt);
		# Assume only one file will have CHR in its name, so
		# we can exit here.
		last if ($pOpt->{c});
	    }
	}
    }
}

# Set up %DNA_BITS and @DNA_SYM for generalizing sets of symbols.

sub setDNABitsAndSymbols {
    %DNA_BITS = qw(A  1 C  2 G  4 T  8
		   M  3 R  5 W  9 
		   S  6 Y 10 K 12
		   B 14 D 13 H 11 V  7
		   N 15 - 16 = 16
		   a 17 c 18 g 20 t 24
		   m 19 r 21 w 25 
		   s 22 y 26 k 28
		   b 30 d 29 h 27 v 23
		   n 31);
    @DNA_SYM = ();
    while (my ($k,$v) = each %DNA_BITS) {
	# Make sure it's been "used as a number" so 
	# bitwise "|" works numerically!
	$DNA_BITS{$k} += 0;
	$DNA_SYM[$v] = $k;
    }
}	      
    
# $g (-o switch) is the file to write to; if undef, we want STDOUT.
# $turnLL (-t switch) is whether to turn the output by 90 degrees.
#   Value is number of nucleotides per line.
# $turnPoly (-P switch): should only be set if $turn is true.

sub openToWrite {
    my ($g,$turnLL,$turnPoly,$tac) = @_;
    # Turning script: will only be used if $turn. It's this
    # script, but with -T 
    my $turner = "$0 -T $turnLL";
    $turner .= " -P" if ($turnPoly);
    my $tacCmd = $tac ? "| tac " : "";
    if (!$g) {
	# No output file specified: print to STDOUT, but pipe via the
	# turning script if -t was given.
	if ($turnLL) {
	    return new FileHandle("$tacCmd| $turner");
	} elsif ($tac) {
	    return new FileHandle($tacCmd);
	} else {
	    return *STDOUT;
	}
    } elsif ($g =~ /\.gz$/) {
	# $g was specified, and is a gzip file. Pipe via turning script
	# if required, then gzip to the file.
	return new FileHandle(($turnLL ? "$tacCmd| $turner " : 0) . "$tacCmd| gzip -c > $g");
    } elsif ($turnLL) {
	# $g was specified, and is not a gzip file. Go via turner.
	return new FileHandle("$tacCmd| $turner > $g");
    } elsif ($tac) {
	return new FileHandle("$tacCmd > $g");
    } else {
	# $g was specified; not gzip; no turner. Write direct to the file.
	return new FileHandle($g,'w');
    }
}

sub alicatRegion {
    my ($pOpt,@files) = @_;
    # Look in the GFF file that should have been specified in the value
    # of -r, and find chromosome, low and high offsets, direction and
    # the identifier also specified in -r. 
    my ($chr,$lo,$hi,$fc,$id) = locateRegion($pOpt->{r});
    # We're going to call alicat1 which expects $pOpt to be populated,
    # so we set the -c, -l and -h switches. We add to -l and -h because
    # they may have been set to non-zero values already, e.g. "-l -100"
    # means "start 100bp before the start of the specified region".
    $pOpt->{c} = $chr;
    $pOpt->{l} += $lo;
    $pOpt->{h} += $hi;
    $pOpt->{R} = $fc eq 'C' ? 1 : undef;
    my $gh = openToWrite($pOpt->{o},$pOpt->{t},$pOpt->{P},$pOpt->{R});
    # Find the (single) file containing the region we want, and call
    # alicat1.
    foreach my $f (@files) {
	if (hasNameComponent($f,$pOpt->{c})) {
	    $pOpt->{c} = "$pOpt->{c}.$id"; # to appear in header
	    alicat1($f,$gh,$pOpt);
	    return;
	}
    }
    # No file for the relevant chromosome appeared on the command line.
    die("Cannot find any file for $pOpt->{c} in @files\n");
}
    
# $pair is the value of the -r switch, of the form id:file, e.g.
#       CUP1:/home/foo/genome.gff
# "id" should be the value of some equation in the last field of
# a GFF line; e.g. a gene name or alias.

sub locateRegion {
    my ($pair) = @_;
    my ($id,$file) = split(/:/,$pair,2);
    # Check the value was in the right format and the file exists.
    if (!defined $file) {
	die("Usage: $0 ... -r identifier:file.gff",1);
    } elsif (! -s $file) {
	die("Missing or empty: $file",1);
    }
    my $fh = openToRead($file);
    while (<$fh>) {
	chomp;
	my @a = split;
	# Equations are separated by ";"
	foreach my $eqn (split(/;/,$a[-1])) {
	    # Each equation should be of the form
	    #   key=val1,val2,...,valN
	    my ($key,$valList) = split(/=/,$eqn,2);
	    my @vals = split(/,/,$valList);
	    # Find a val that equals $id.
	    foreach my $val (@vals) {
		if ($val eq $id) {
		    return ($a[0],$a[3],$a[4],$a[6] eq '-' ? 'C' : 'F',$id);
		}
	    }
	}
    }
    die("Cannot find any line matching $id in $file",1);
}

# The file may be unspecified (meaning STDIN), or if specified
# may be gzipped or not.

sub openToRead {
    my ($f) = @_;
    if (!$f) {
	return *STDIN;
    } elsif ($f =~ /\.gz$/) {
	return new FileHandle("gunzip -c $f |");
    } else {
	return new FileHandle($f,'r');
    }
}

# Return true if any component of the pathname $f equals $id.

sub hasNameComponent {
    my ($f,$id) = @_;
    foreach my $comp (split(/\//,$f)) {
	return 1 if ($comp eq $id || $comp =~ /^$id\./);
    }
    return 0;
}

# "accept" values in alicat1 are 1 or 0 for consensus case.  For
# non-consensus, "accept" value for each position (0 being cols 3 and
# 4) is the position we map to, and acceptBack goes the other way.

# $f: file to read from
# $gh: output file handle
# $pOpt: command-line options (which may have been adjusted by alicatRegion).

sub alicat1 {

    my ($f,$gh,$pOpt) = @_;

    # First set up values from switches ...
    # Look up chromosome name, low and high offset.
    my ($chr,$lo,$hi) = ($pOpt->{c},$pOpt->{l},$pOpt->{h});
    # Always print, even if just gaps
    my $always = $pOpt->{A};
    # $polyOnly: if true, print only polymorphic positions.
    # $qTh: quality threshold for counting a polymorphism and for
    # not printing a nucleotide as "N".
    my ($polyOnly,$qTh) = ($pOpt->{p},$pOpt->{q});
    # $sep: if true, print nucleotides in one "|"-separated block
    # and qualities in another; otherwise, print interleaved, as in
    # the input files.
    my $sep = $pOpt->{s};
    # Translate numeric gap symbols to "-" UNLESS -g was set.
    my $trans = !$pOpt->{g};
    # Patterns for including/excluding strain names.
    my ($incPats,$excPats) = ($pOpt->{i},$pOpt->{x});
    # Whether to turn strain names vertical. We do this
    # if either -v or -t was specified.
    my $vertHdr = $pOpt->{v} || $pOpt->{t};
    # "-a" means use alignment coordinates, i.e. increment at all
    # alignment lines, even when ref is a gap there.
    my $incrementOffsetAtGaps = $pOpt->{a};
    # $tensQualities is true if qualities should be printed as
    # numbers: 0 to 9 as 0, 10 to 19 as 1, etc.
    my $tensQualities = $pOpt->{y};
    # Use dots and colons for repeated values (see -d)
    my $dotsAndColons = $pOpt->{d};
    # Reverse-complement nucleotides.
    my $revcomp = $pOpt->{R};

    # Offset: start at 0 if -z, else 1 (this value gets
    # incremented before anything happens). Increment at reference
    # nucleotides only (no -a) or at all data lines (-a).
    my $n = $pOpt->{z} ? -1 : 0;
    # Delta-offset: increment for each reference gap, unless -a
    my $d=0;

    # Whether protein sequence: look for line starting ">... gene ..."
    my $isProtein;

    # Not documented: called from imputed2html in web site generation
    # (which needs work). Better would be: if -r was specified, and it's
    # for a CDS, use that information. And allow -r to _just_ be the GFF
    # file, then we'll print everything.
    my $pCodons;
    if ($pOpt->{G} && !$pOpt->{a}) {
	my $ch = openToRead($pOpt->{G});
	while (<$ch>) {
	    chomp;
	    my ($off,@rest) = split;
	    $off -= 2 if ($rest[0] eq '-'); # shift to LO end of codon even when reverse
	    push @{$pCodons->{$off}},\@rest;
	}
    }

    # In consensus (sorted/imputed) files, $accept[$i] is 1 for strains we want,
    # else 0, and @acceptBack is not used. In non-consensus (layout/rqpw),
    # $accept[$i] is the output position for input position $i, or undef if
    # position is not wanted, and @acceptBack is the inverse of @accept.
    my (@accept,@acceptBack);
    # $expectedLineLength is the number of characters in a line when all strains
    # have data -- only defined for consensus files.
    # $wanted{$str} is whether we want strain $str. $anno[$i] is the pending annotation
    # for INPUT position $i, filled when we're still in quiet mode; we'll print all pending
    # annotations when we come out of quiet mode.
    my ($expectedLineLength,%wanted,@anno);
    # If $isConsensus, it's a consensus (per strain) file, rather than a 
    # per-read file.
    my $isConsensus = !($OPT{L} || $f =~ /layout/ || $f =~ /rqpw/);

    # Let's get reading ...
    my $fh = openToRead($f);
    
    if ($isConsensus) {
	# For a consensus file, the first line should give the strain names.
	chomp (my $hdr = <$fh>);
	# All the strains in the header (plus REF at the start):
	my @strains = ('REF',split(" ",substr($hdr,1)));
	# $accept[$i] is 1 if we should print $strains[$i]
	@accept = applyIncAndExc($incPats,$excPats,\@strains);
	# Expected input line length:
	$expectedLineLength = scalar(@strains)*2;
	# Print the header, perhaps vertically.
	printConsensusHeader($vertHdr,\@strains,\@accept,$sep,$chr,$gh);
    }

    # Now loop round, reading every line.
    while (<$fh>) {
	last if (defined $hi && $n >= $hi); # gone far enough
	chomp;
	# We're "quiet" if we haven't reached offset $lo yet.
	my $quiet = defined $lo && $n < $lo-1;
	# Lengthen the line to the expected length if necessary
	# (for consensus files only).
	if ($isConsensus) {
	    my $dl = $expectedLineLength-length($_);
	    $_ .= ' ' x $dl if ($dl > 0);
	}
	if (/^>/) {
	    if ($isConsensus) {
		if (/\sgene\s/) {
		    $isProtein = 1;
		    $n = $d = 0;
		}
		# Annotation line. Shouldn't (often) happen for a consensus file,
		# as we'll already have dealt with the header, but anyway we print
		# the line if we're printing.
		printConsensusAnnoLine($chr,$_,$gh) unless ($quiet || $vertHdr);
	    } else {
		# For non-consensus, we either print, or we save it up for possible
		# printing later.
		printNonConsensusAnnoLine($_,$incPats,$excPats,\%wanted,
					  \@accept,\@acceptBack,\@anno,$quiet,$gh);
	    }
	} elsif ($quiet) {
	    # Data line, but we're quiet. Just increase $n and $d (don't hand over $gh)
	    ($n,$d) = printData($_,undef,$polyOnly,$incrementOffsetAtGaps,$n,$d,undef);
	    # If we're just about to start printing data, i.e. $quiet will be false
	    # next time, then print all pending annotations.
	    if ($n == $lo-1) {
		dischargeAnnotations(\@accept,\@acceptBack,\@anno,$gh) 
	    }
	    # Relinquish any newly-space columns, but don't bother to build the output line
	    # as we won't print it.
	    contractNonConsensusDataLine(\@accept,\@acceptBack,$_,0) unless ($isConsensus);
	} else {
	    # We're going to print ...
	    if ($isConsensus) {
		# Just keep the pairs from the "accept" positions ...
		$_ = contractConsensusDataLine(\@accept,$_);
	    } else {
		$_ = contractNonConsensusDataLine(\@accept,\@acceptBack,$_,1);
	    }
	    # We do this before finding the consensus!
	    $_ = reverseComplementNucleotides($_) if ($revcomp);
	    # Find the consensus value (generalization) of all the nucleotides we're going to print.
            # It'll be bounded by spaces if it's a simple value, or "[ ]" if it's a generalization.
	    my $cons = consensusValue($_,$qTh,$isProtein);
	    # If consensus is a gap, we don't print the line. This 
	    # can happen if we have exclusions and only the excluded strains have nucleotides.
	    unless ((!$always) && ($cons eq ' - ' || $cons eq ' = ')) {
		# Apply various transformations as defined above.
		$_ = translateDigitGaps($_) if ($trans);
		$_ = introduceDotsAndColons($_) if ($dotsAndColons);
		$_ = poorNucleotidesToN($_,$qTh) if ($qTh > 0);
		$_ = applyTensQualities($_) if ($tensQualities);
		$_ = separateNucsAndQuals($_) if ($sep);
		# And print the data for the transformed line.
		($n,$d) = printData($_,$cons,$polyOnly,$incrementOffsetAtGaps,$n,$d,$gh,$pCodons);
	    }
	}
    }
}

sub poorNucleotidesToN {
    my ($line,$th) = @_;
    my $thc = $th >= 93 ? '~' : chr($th+33);
    my $iLim = length($line);
    for (my $i=0; $i<$iLim; $i+=2) {
	my $v = substr($_,$i,1);
	if ($v =~ /[ACGTacgt]/) {
	    my $qc = substr($_,$i+1,1);
	    if ($qc ne ' ' && $qc lt $thc) {
		substr($_,$i,1) = ($v =~ /[ACGT]/ ? 'N' : 'n');
	    }
	}
    }
    return $_;
}	

# Given (comma-separated) $incPats and $excPats, return list of 1's and 0's for
# whether to keep each strain in @$pStrains.

sub applyIncAndExc {
    my ($incPats,$excPats,$pStrains) = @_;
    my @accept;
    foreach my $str (@$pStrains) {
	push @accept,wantThisStrain($incPats,$excPats,$str);
    }
    return @accept;
}

sub wantThisStrain {
    my ($incPats,$excPats,$str) = @_;
    my @incPats = split(',',$incPats);
    my @excPats = split(',',$excPats);

    return 1 if ($str eq 'REF'); # Always want REF.

    # Step 1: $got is 1 if there are no @incPats or if the strain is
    # REF. Otherwise, $got is 1 if the whole of the strain name
    # matches one of the patterns.
    my $got = 0;
    if (@incPats) {
	foreach my $pat (@incPats) {
	    if ($str eq $pat || $str =~ /^$pat$/) {
		$got = 1;
		last;
	    }
	}
    } else {
	$got = 1;
    }
    return 0 unless ($got); # $str is not included.
    # Now see if it's to be excluded.
    foreach my $pat (@excPats) {
	if ($str eq $pat || $str =~ /^$pat$/) {
	    return 0; # yes,it is.
	}
    }
    return 1; # Still OK: included and not excluded.
}

# In vertical-header mode, print a consensus-file header like 
sub printConsensusHeader {
    my ($vertHdr,$pStrains,$pAccept,$sep,$chr,$gh) = @_;
    if ($vertHdr) {
	printVerticalHeader($pStrains,$pAccept,$sep,$gh);
    } else {
	printHorizontalHeader($pStrains,$pAccept,$chr,$gh);
    }
}

# Print ">", then the chr if any, then the accepted strains.

sub printHorizontalHeader {
    my ($pStrains,$pAccept,$chr,$gh) = @_;
    print $gh ">";
    print $gh "$chr " if ($chr);
    my $sepr = '';
    for (my $k=0; $k<=$#$pStrains; $k++) {
	if ($pAccept->[$k]) {
	    print $gh "$sepr$pStrains->[$k]";
	    $sepr = ' ';
	}
    }
    print $gh "\n";
}

sub printVerticalHeader {
    my ($pStrains,$pAccept,$sep,$gh) = @_;
    my $k=0; # column number
    my $cons = 'Cons';
    while (1) {
	# 13 spaces, then a character of "Cons" if there any left, then 3 spaces.
	my $line = sprintf("%13s%1s%3s","",$k >= 4 ? ' ' : substr($cons,$k,1),"");
	# Print the k'th character of every accepted strain name (or ' ' if shorter than $k).
	# If not $sep, we'll have two columns per strain (for nuc and qual).
	$line .= printStrainRow($k,$pStrains,$pAccept,$sep);
	# If $sep, i.e. -s switch, print 3 spaces, then do it all again.
	if ($sep) {
	    $line .= '   ';
	    $line .= printStrainRow($k,$pStrains,$pAccept,$sep);
	}
	# May as well get rid of any trailing spaces.
	$line =~ s/  +$//;
	print $gh "$line\n";
	last unless ($line);
	$k++;
    }
}

sub printStrainRow {
    my ($k,$pStrains,$pAccept,$sep) = @_;
    my $line;
    my @acc = @$pAccept;
    foreach my $str (@$pStrains) {
	next unless (shift @acc);
	if (length($str) > $k) {
	    $line .= substr($str,$k,1);
	} else {
	    $line .= ' ';
	}
	$line .= ' ' unless ($sep);
    }
    return $line;
}

# Print the line if either we have a $chr to print or the line has some content.
# This shouldn't usually get called, though.

sub printConsensusAnnoLine {
    my ($chr,$line,$gh) = @_;
    if ($chr) {
	# Smuggle in dir name
	printf($gh ">%s %s\n",$chr,substr($line,1));
    } elsif ($line =~ /^>.*\S/) {
	# Print if we have something other than a space after the ">"
	print $gh "$line\n";
    }
}

# In a non-consensus (e.g. layout) file, we first map the position (where
# column = 2*(position+1)) onto a new value appropriate for the strains we're
# interested in, assuming we do want this strain. Then we either print the
# resulting trio or store it in $pAnno->[$pos] in case it's still in force
# when we leave "quiet" mode.

sub printNonConsensusAnnoLine {
    my ($line,$incPats,$excPats,$pWanted,$pAccept,$pAcceptBack,$pAnno,$quiet,$gh) = @_;
    my ($name,$pos,$iter) = split(" ",substr($line,1));
    if ($name eq 'REF') {
	# do nothing for now, as these numbers are often wrong!
	#print $gh "$line\n" unless ($quiet);
    } else {
	my ($str) = split(/-/,$name);
	if (wantStrain($str,$incPats,$excPats,$pWanted)) {
	    # $tgt is position to use in the output. Store it or print it.
	    my $tgt = introduceRead($pos,$pAccept,$pAcceptBack);
	    my $trio = "$name $tgt $iter";
	    if ($quiet) {
		$pAnno->[$pos] = $trio;
	    } else {
		print $gh ">$trio\n";
	    }
	}
    }
}

sub introduceRead {
    my ($pos,$pAccept,$pAcceptBack) = @_;
    # If something was already in this position, wipe it out.
    if (defined $pAccept->[$pos]) {
	$pAcceptBack->[$pAccept->[$pos]] = undef;
	$pAccept->[$pos] = undef;
    }
    # Find first unoccupied position and use it.
    for (my $i=0; 1; $i++) {
	unless (defined $pAcceptBack->[$i]) {
	    $pAccept->[$pos] = $i;
	    $pAcceptBack->[$i] = $pos;
	    return $i;
	}
    }
}

sub wantStrain {
    my ($str,$incPats,$excPats,$pWant) = @_;
    my $v = $pWant->{$str};
    unless (defined $v) {
	$v = $pWant->{$str} = wantThisStrain($incPats,$excPats,$str);
    }
    return $v;
}

sub printData {
    my ($line,$cons,$polyOnly,$incrementOffsetAtGaps,$n,$d,$gh,$pCodons) = @_;

    # Polymorphic line if the consensus symbol is bracketed.
    my $isPoly = substr($cons,0,1) eq '[';

    # We want to print the line if we have an output stream and we either want
    # all lines or this line is polymorphic.
    my $active = $gh && !($polyOnly && !$isPoly);

    # Increment the offset (either $n or $d), and print the offset, consensus and data line.
    if ((!$incrementOffsetAtGaps) && $line =~ /^[-=_0-9a-z ]/) { # gap or ref-genome extension
	$d ++;
	printf($gh "%7d.%s %s |%s|",$n,makeTrio($d),$cons,$line) if ($active);
    } else {
	$d=0;
	$n++;
	printf($gh "%7d     %s |%s|",$n,$cons,$line) if ($active);
    }

    # Print the mutation effect if we have codons available.
    if ($active) {
	if ($pCodons) {
	    print $gh " " . mutationEffect($isPoly,substr($cons,1,1),$pCodons,$n);
	    if ($d == 0) {
		if ($pCodons->{$n}) {
		    my $sep = ' ';
		    foreach my $pA (@{$pCodons->{$n}}) {
			printf($gh "%s%s",$sep,join(" ",@$pA));
			$sep = '; ';
		    }
		} elsif ($pCodons->{$n-3}) {
		    printf($gh " %s Intergenic",$pCodons->{$n-3}[0][0]);
		}
	    }
	}
	print $gh "\n";
    }
    return ($n,$d);
}

# Effect values:
#   'X'  external to coding sequence
#   'S'  synonymous (non-coding)
#   'N'  non-synonmyous
#   'P'  phase change (indel)
#   'T'  stop introduction/removal

sub mutationEffect {
    my ($isPoly,$consChar,$pCodons,$n) = @_;
    return ' ' unless ($isPoly);
    return 'P' if ($consChar =~ /[a-z]/);
    my $maxEffect = 0;
    foreach my $m ($n-2 .. $n) {
	if (defined $pCodons->{$m}) {
	    foreach my $pA (@{$pCodons->{$m}}) {
		$maxEffect =  max($maxEffect,mutationEffect1($n-$m,$pA,$consChar));
	    }
	}
    }
    my @eff = ('X','S','N','T');
    return $eff[$maxEffect];
}
   
sub max {
    my ($x,$y) = @_;
    if (defined $y && $y > $x) {
	return $y;
    } elsif (defined $x) {
	return $x;
    } else {

	return $y;
    }
}

sub reverseComplementNucleotides {
    my ($line) = @_;
    my $iLim = length($line);
    for (my $i=0; $i<$iLim; $i+=2) {
	substr($line,$i,1) =~ tr/ACGTacgt/TGCAtgca/;
    }
    return $line;
}

sub mutationEffect1 {
    my ($pos,$pA,$cons) = @_;
    my ($dir,$codon) = @$pA;
    $pos = 2-$pos if ($dir eq '-');
    my $aa = $TRANSLATION{$codon};
    my $bits = $DNA_BITS{$cons};
    my $maxEffect = 0;
    foreach my $val (qw(A C G T)) {
	if (($DNA_BITS{$val} & $bits) && substr($codon,$pos,1) ne $val) {
	    my $codon2 = $codon;
	    my $val2 = $val;
	    $val2 =~ tr/ACGT/TGCA/ if ($dir eq '-');
	    substr($codon2,$pos,1) = $val2;
	    my $aa2 = $TRANSLATION{$codon2};
	    my $effect;
	    if ($aa2 eq $aa) {
		$effect = 1;
	    } elsif (($aa2 eq 'STP') eq ($aa eq 'STP')) {
		$effect = 2;
	    } else {
		$effect = 3;
	    }
	    $maxEffect = max($maxEffect,$effect);
	}
    }
    return $maxEffect;
}
	

sub dischargeAnnotations {
    my ($pAccept,$pBack,$pAnno,$gh) = @_;
    foreach my $i (@$pBack) {
	if (defined $i && defined $pAccept->[$i]) {
	    print $gh ">$pAnno->[$i]\n";
	}
    }
}

sub makeTrio {
    my ($d) = @_;
    my $dd;
    if ($d < 1000) {
	$dd = sprintf("%3d",$d);
    } else {
	my $h = int($d/100)-10;
	my $u = $d%100;
	my $hc;
	if ($h < 26) {
	    $hc = chr(65+$h);
	} elsif ($h < 52) {
	    $hc = chr(71+$h);
	} else {
	    $hc = '*';
	}
	$dd = sprintf("%s%2d",$hc,$u);
    }
    $dd =~ tr/ /0/;
    return $dd;
}

# Given "accept" and "acceptBack" arrays, and an input data line,
# return the contracted (output) line.

sub contractNonConsensusDataLine {
    my ($pAcc,$pBack,$line,$really) = @_;
    # Output line; start with reference pair.
    my $ans = substr($line,0,2);
    # Length of output line.
    my $lenAns = 2;
    # Current input column for position $pos (0 -> 2, 1 -> 4 etc).
    my $col = 0;
    my $len = length($line);
    for (my $pos=0; $col<=$len; $pos++) {
	$col += 2;
	if (substr($line,$col,1) eq ' ') {
	    # If we've got a space at this position, relinquish it if
	    # it's currently assigned.
	    if (defined $pAcc->[$pos]) {
		$pBack->[$pAcc->[$pos]] = undef;
		$pAcc->[$pos] = undef;
	    }
	} elsif ($really && defined $pAcc->[$pos]) {
	    # Extend the output line: $dst is the destination column.
	    my $dst = 2*(1+$pAcc->[$pos]);
	    # Extend with spaces if needed.
	    if ($lenAns <= $dst) {
		$ans .= ' ' x ($dst+2-$lenAns);
		$lenAns = $dst+2;
	    }
	    # Fill the column (and the one after it, for the quality).
	    substr($ans,$dst,2) = substr($line,$col,2);
	}
    }
    # Relinquish positions beyond current line length
    for (my $pos=$len/2; $pos <= $#$pAcc; $pos++) {
	if (defined $pAcc->[$pos]) {
	    $pBack->[$pAcc->[$pos]] = undef;
	    $pAcc->[$pos] = undef;
	}
    }
    # Shorten the @accept array by popping any undef values.
    while (@$pAcc && !defined $pAcc->[-1]) {
	pop @$pAcc;
    }
    # Return the output line (should only be used if $really).
    return $ans;
}
    
sub contractConsensusDataLine {
    my ($pAccept,$line) = @_;
    my $ans;
    foreach my $i (0 .. $#$pAccept) {
	if ($pAccept->[$i]) {
	    $ans .= substr($line,$i*2,2);
	}
    }
    return $ans;
}

sub consensusValue {
    my ($line,$qTh,$isProtein) = @_;
    my $kMax = length($line)-2;
    my $cons = substr($line,0,1);
    $cons =~ tr/_a-z/-A-Z/;
    $cons =~ s/[=0-9]/-/g;	    
    my $cons1 = $cons;
    my $qc;
    $qc = chr($qTh+33) if (defined $qTh);
    for (my $k=2; $k<=$kMax; $k+=2) {
	if (defined $qc) {
	    my $qk = substr($line,$k+1,1);
	    next if ($qk ne ' ' && $qk lt $qc);
	}
	$cons = consensusValue1($cons,substr($line,$k,1),$isProtein);
	#printf("'%s' '%s'\n",substr($line,$k,1),$cons);
    }
    return $cons1 eq $cons ? " $cons " : "[$cons]";
}

sub consensusValue1 {
    my ($v1,$v2,$isProtein) = @_;
    if ($v2 eq ' ' || $v1 eq $v2 || $v1 eq '*') {
	return $v1;
    } elsif ($isProtein) {
	if ($v2 eq '-' || $v2 eq '+') {
	    $v1 =~ tr/A-Z/a-z/;
	    return $v1;
	} elsif ($v1 eq '-' || $v1 eq '+') {
	    $v2 =~ tr/A-Z/a-z/;
	    return $v2;
	} else {
	    my $u1 = $v1;
	    $u1 =~ tr/_a-z/-A-Z/;
	    $u1 =~ s/[0-9]/-/g;	   
	    my $u2 = $v2;
	    $u2 =~ tr/_a-z/-A-Z/;
	    $u2 =~ s/[0-9]/-/g;	   
	    # Preserve lower case:
	    return ($u2 eq $u1) ? ($u1 eq $v1 ? $v2 : $v1) : '*';
	}
    } else {
	$v2 =~ tr/_a-z/-A-Z/;
	$v2 =~ s/[0-9]/-/g;	   
	return $DNA_SYM[$DNA_BITS{$v1} | $DNA_BITS{$v2}];
    }
}

sub translateDigitGaps {
    my ($line) = @_;
    my $kMax = length($line)-2;
    for (my $k=0; $k<=$kMax; $k+= 2) {
	substr($line,$k,1) =~ s/[0-9]/-/;
    }
    return $line;
}

sub introduceDotsAndColons {
    my ($line) = @_;
    my $rColon = substr($line,0,1);
    (my $rDot = $rColon) =~ tr/-A-Z/_a-z/;
    my $iMax = length($line)-2;
    for (my $i=2; $i<=$iMax; $i+=2) {
	my $v = substr($line,$i,1);
	if ($v eq $rColon) {
	    substr($line,$i,1) = ':';
	} elsif ($v eq $rDot) {
	    substr($line,$i,1) = '.';
	}
    }
    return $line;
}

sub applyTensQualities {
    my ($line) = @_;
    my $iMax = length($line)-1;
    for (my $i=3; $i<=$iMax; $i+=2) {
	my $qc = substr($line,$i,1);
	unless ($qc eq ' ') {
	    my $q = ord($qc) - 33;
	    substr($line,$i,1) = int($q/10);
	}
    }
    return $line;
}


sub separateNucsAndQuals {
    my ($line) = @_;
    my $kMax = length($line)-2;
    my ($nucs,$quals);
    for (my $k=0; $k<=$kMax; $k+= 2) {
	$nucs .= substr($line,$k,1);
	$quals .= substr($line,$k+1,1);
    }
    return "$nucs| |$quals";
}		

# Input: output of alicat.pl

sub aliturn {
    my ($lineLen,$polyOnly) = @_;
    my @pfx;
    my @lines;
    while (<>) {
	chomp;
	s/ +$//;
	my $kMax = length($_)-1;
	if ($kMax < 0 && !@pfx) {
	    @pfx = @lines;
	    @lines = ();
	}
	for (my $k=0; $k<=$kMax; $k++) {
	    if ($k > 0 && !$lines[$k]) {
		$lines[$k] = ' ' x length($lines[0]);
	    }
	    $lines[$k] .= substr($_,$k,1);
	}
	for (my $k=$kMax+1; $k<=$#lines; $k++) {
	    $lines[$k] .= ' ';
	}
	if (length($lines[0]) >= $lineLen) {
	    printTurnedLines($lineLen,\@pfx,\@lines,$polyOnly);
	    @lines = ();
	}
    }
    printTurnedLines($lineLen,\@pfx,\@lines,$polyOnly);
}

sub printTurnedLines {
    my ($lineLen,$pPfx,$pLines,$polyOnly) = @_;
    # Put spaces into line number lines
    foreach my $k (0 .. $#$pLines) {
	if ($pLines->[$k] =~ /^ *[0-9]+$/) {
	    for (my $i=length($pLines->[$k])-1; $i>=1; $i--) {
		if (substr($pLines->[$k],$i,1) eq substr($pLines->[$k],$i-1,1)) {
		    substr($pLines->[$k],$i,1) = ' ';
		}
	    }
	} elsif ($pLines->[$k] =~ /\S/) {
	    last;
	}
    }
    my $lastBlank = 1;
    print ' ' x (2+length($pPfx->[0]));
    print '-' x $lineLen;
    print "\n";
    my $postRef = 0;
    my $skipQual;
    foreach my $k (0 .. $#$pLines) {
	if ($skipQual) {
	    $skipQual = 0;
	    #print "## SKIPQ '$pPfx->[$k]'\n";
	    next;
	}
	$skipQual = 0;
	$postRef = 1 if ($pPfx->[$k] =~ / REF /);
	if ($pPfx->[$k] =~ /\S/) {
	    if ($polyOnly && $postRef) {
		if ($pLines->[$k] !~ /[-=_0-9A-Za-z]/) {
		    #print "## SKIP '$pPfx->[$k]'\n";
		    $skipQual = 1;
		    next;
		}
	    }
	} else {
	    $pPfx->[$k] = ' ' x length($pPfx->[1]);
	}
	my $line = "$pPfx->[$k]$pLines->[$k]";
	$line =~ s/ +$//;
	if ($line =~ /[^ |]/) {
	    print "$line\n";
	    $lastBlank = 0;
	} else {
	    print "\n" unless ($lastBlank);
	    $lastBlank = 1;
	}
    }
}
