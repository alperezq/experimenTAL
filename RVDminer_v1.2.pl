#!/usr/bin/perl
use strict;
use warnings;
use Bio::Perl;
use Bio::DB::Fasta;
use Bio::SeqIO;
use List::MoreUtils qw(uniq);
use List::Util qw( min max );
#use String::Approx 'aindex';

#work on getting strain and species info, and on better matching the strings

my $help = '
usage: RVDminer.pl Input Output Genebankfile(optional)
Input = fasta file with nucleotide sequences of TAL effectors
';


my %db;
my $obj;

#create formatted fasta database object
if ( $ARGV[0] ) {
	$obj = tie %db, 'Bio::DB::Fasta', "$ARGV[0]";
}

else { print "$help"; exit; }


#Get information on strains and species or the way to get this information 

my $genebank;
my $single_specie;
my $single_strain;

my $multiple_species = "F";
my $multiple_strains = "F";
my $delim1;my $field1;
my $delim2;my $field2;


print "Analizing multiple (M) or single species (S)?\n";
my $in2 = <STDIN>; chomp $in2;
if ($in2 eq "S") {print "insert species name:\n"; $single_specie = <STDIN>; chomp $single_specie}
elsif ($in2 eq "M") 				{
print "If species are included in the sequence header, state delimiter and field: example: | 3 \n";
my $in3 = <STDIN>; chomp $in3; $delim1 = "\\".(split(/\s/,$in3))[0]; $field1 = (split(/\s/,$in3))[1]; $multiple_species = "T";
						}
print "Analizing multiple (M) or single strains (S)?\n";
my $in4 = <STDIN>; chomp $in4;
if ($in4 eq "S") {print "insert strain name:\n"; $single_strain = <STDIN>; chomp $single_strain}
elsif ($in4 eq "M") 				{
print "If strainss are included in the sequence header, state delimiter and field: example: | 3 \n";
my $in5 = <STDIN>; chomp $in5; $delim2 = "\\".(split(/\s/,$in5))[0]; $field2 = (split(/\s/,$in5))[1]; $multiple_strains = "T";
						}


#read features from genebank file to be included in annotation (optional)

my $seqio_object;
my %strains;
my %species;

#create output files

open my $OUT1, "> $ARGV[1]"."_Nterminal_nuc.FASTA";
open my $OUT2, "> $ARGV[1]"."_Cterminal_nuc.FASTA";
open my $OUT3, "> $ARGV[1]"."_Repeats_nuc.FASTA";
open my $OUT4, "> $ARGV[1]"."_RVDs_nuc.csv";
print $OUT4 "ID\tRVDseq\tSpecie\tStrain\tAccesion\tgi\tComments\n";
#open my $OUT5, "> $ARGV[1]"."_0Repeat_nuc.FASTA";
open my $OUT7, "> $ARGV[1]"."_TALS_nuc.FASTA";
open my $OUT6, "> $ARGV[1]"."_TALS_aa.FASTA";
open my $OUT8, "> $ARGV[1]"."_Repeats_aa.FASTA";
open my $OUT9, "> $ARGV[1]"."_TALs.gff3";

#####array of possible starting nucleotides in a TAL repetition

my @Ini_repeat;
open my $IN1, "<Initial_TAL_strings.txt";

while ( my $line = <$IN1> ) {
	chomp $line;
	push @Ini_repeat, $line;
}

close $IN1;

while ( my ( $id, $seq ) = each %db )    #explore each sequence in the database

{

	#get information from IDs
	my $header = $obj->header($id);
	print "$header\n";
	#print "$id\n";
	my $seqaccession = "$header";
	my $seqgi = "-";
	my $aaseq = translate_as_string ($seq);
	my $aaseq2 = translate_as_string (substr $seq, 1);
	my $aaseq3 = translate_as_string (substr $seq, 2);
	my $revcom = reverse_complement_as_string($seq);
	my $aaseq4 = translate_as_string ($revcom);
	my $aaseq5 = translate_as_string (substr $revcom, 1);
	my $aaseq6 = translate_as_string (substr $revcom, 2);
	my $seqlen = length ($seq);

	#$seqaccession = (split (/\,/, (split(/ /,$id))[1]))[0];
	if ($header =~ "[\=]") {
	$seqaccession = (split (/[\=.]/, (split(/\[/,$header))[1]))[1];
	$seqaccession =~ s/] //g;
	$header = $seqaccession;
	#$seqaccession = (split(/\>\|/,$header))[1];
	#print "SEQACCS $seqaccession\n";
	}
		
	if ($multiple_species eq "T") {$species{$seqaccession} = (split (/$delim1/, $header))[$field1];}
	if ($multiple_strains eq "T") {$strains{$seqaccession} = (split (/$delim2/, $header))[$field2];}
	

	#my $sxs = (split(/\|/,$header))[0]; 
	#$header = "$sxs";

	
	my @position_aa_array = ($seqlen);
	my @minus_position_aa_array = ($seqlen);


	#look for each possible TAL repetition initiation in all six frames

	foreach (@Ini_repeat) {
		my $ini_motif = $_;

		my $aaindexresult = 0;
		my $aaindexresult2 = 0;
		my $aaindexresult3 = 0;
		my $aaindexresult4 = 0;
		my $aaindexresult5 = 0;
		my $aaindexresult6 = 0;
		my $offset = 0;


		while ( $aaindexresult != -1 ) {
			$aaindexresult = index($aaseq, $ini_motif, $offset);		
			$offset = $aaindexresult + 1;		
			if ( $aaindexresult > -1) { push @position_aa_array, ($aaindexresult*3); }}
		while ( $aaindexresult2 != -1 ) {
			$aaindexresult2 = index($aaseq2, $ini_motif, $offset);		
			$offset = $aaindexresult2 + 1;		
			if ( $aaindexresult2 > -1) { push @position_aa_array, ($aaindexresult2*3+1); }}
		while ( $aaindexresult3 != -1 ) {
			$aaindexresult3 = index($aaseq3, $ini_motif, $offset);		
			$offset = $aaindexresult3 + 1;		
			if ( $aaindexresult3 > -1) { push @position_aa_array, ($aaindexresult3*3+2); }}
	#do it for the oposite strand
		while ( $aaindexresult4 != -1 ) {
			$aaindexresult4 = index($aaseq4, $ini_motif, $offset);		
			$offset = $aaindexresult4 + 1;		
			if ( $aaindexresult4 > -1) { push @minus_position_aa_array, ($aaindexresult4*3);  }}
		while ( $aaindexresult5 != -1 ) {
			$aaindexresult5 = index($aaseq5, $ini_motif, $offset);		
			$offset = $aaindexresult5 + 1;		
			if ( $aaindexresult5 > -1) { push @minus_position_aa_array, ($aaindexresult5*3+1); }}
		while ( $aaindexresult6 != -1 ) {
			$aaindexresult6 = index($aaseq6, $ini_motif, $offset);		
			$offset = $aaindexresult6 + 1;		
			if ( $aaindexresult6 > -1) { push @minus_position_aa_array, ($aaindexresult6*3+2); }}
				
		

	
	}

	#create array with uniq positions where the the initial TAL string was found
	my @POS = uniq( sort { $a <=> $b } (@position_aa_array) ); 

	#create array with uniq positions where the the initial TAL string was found in the minus strand
	my @minPOS = uniq( sort { $a <=> $b } (@minus_position_aa_array) ); 

	#Find number of TALs and create arrays for each one
	
	my @AoA; #array of array for positions

	#my @array = @{$AoA[0]};
	#$AoB[$i] = [ @array ];
	
	#create arrays of positions for each TAL

	my $TALi = 0;
	for ( my $i = 1 ; $i <= $#POS ; $i++ ) {
	my $length = $POS[$i] - $POS[ $i - 1 ];
	if ($length < 800) {push @{$AoA[$TALi]}, $POS[ $i - 1 ];}
	else {push @{$AoA[$TALi]}, $POS[ $i - 1 ]; $TALi++; print "FOUnd another TAL!!\n"}
	}
	
	my @AoB;
	#minus strand
	$TALi = 0;
	for ( my $i = 1 ; $i <= $#minPOS ; $i++ ) {
	my $length = $minPOS[$i] - $minPOS[ $i - 1 ];	
	if ($length < 800) {push @{$AoB[$TALi]}, $minPOS[ $i - 1 ]}
	else {push @{$AoB[$TALi]}, $minPOS[ $i - 1 ]; $TALi++; print "FOUnd another TAL!!\n"}
	}

my @Ntermmotifs = ("ATGGATCCC","ATGAGAATA");
my @STOPmot = ("TAG","TAA","TGA");

for (my $j = 0; $j <= $#AoA; $j++) { #FOR EACH TAL in the plus strand
	my @comments; #comments to add about each sequence
	my $TALlen = $#{$AoA[$j]};
	if ($TALlen > 0) {
	print "TALlength = $TALlen\n";
	#extract repetitions and RVDs
	
	my @Repeat_len;
	my @Repeat_array;
	my @AArepeat;
	my @RVD_array;
	my $RVD;
	
	#find repeat length
	for ( my $i = 1 ; $i <= $TALlen ; $i++ ) {
		my $length = $AoA[$j][$i] - $AoA[$j][ $i - 1 ];
		push @Repeat_len, $length;
		
	}
	#find most common repetition length
		my(%count);
		foreach my $value (@Repeat_len) {
    		$count{$value}++;
		}
		my $common_len = (sort {$count{$b} <=> $count{$a}} @Repeat_len)[0];
		my $min_len = (sort {$a <=> $b} @Repeat_len)[0];
		my $max_len = (sort {$a <=> $b} @Repeat_len)[$#Repeat_len];
	

	for ( my $i = 0 ; $i <= $#Repeat_len ; $i++ ) {
		my $repeat = substr ($seq, $AoA[$j][$i], $Repeat_len[$i]);
		push @Repeat_array, $repeat;
		my $aarepeat = translate_as_string ($repeat);
		push @AArepeat, $aarepeat;
		#print "$aarepeat\n";
		if ($Repeat_len[$i] == $common_len) {$RVD = substr ($aarepeat, 11, 2); push @RVD_array, $RVD;}
		#add asterisk to RVD if repeat is one aminoacid shorter
		elsif ($Repeat_len[$i] == $common_len - 3) {$RVD = substr ($aarepeat, 11, 1); push @RVD_array, $RVD."*"} 
		else {$RVD = substr ($aarepeat, 11, 2); push @RVD_array, $RVD; push @comments, "aberrant repeat;"; print "Aberrant or unrecognized repeat : $aarepeat\n"}
						}
	#get last repeat		
		my $last_repeat = substr ($seq, $AoA[$j][$TALlen], $common_len/2);
		push @Repeat_array, $last_repeat;
		my $aalast_repeat = translate_as_string ($last_repeat);
		print "$aalast_repeat\n";
		push @AArepeat, $aalast_repeat;
		push @RVD_array, substr ($aalast_repeat, 11, 2);

	#find N terminal


	my $ATG = -1;
	my $Nregion;
	my $Nseq;

	my $Repstart = $AoA[$j][0];
	if ($Repstart - 1500 > 0) {$Nregion = $Repstart - 1500; $Nseq =  substr $seq, $Nregion, 1500}
		else {$Nregion = 0; $Nseq = substr $seq, 0, $Repstart}
	print "REPstart = $Repstart\n";
	foreach (@Ntermmotifs) {
		my $motif = $_;
		my $indexresult = 0;
		$indexresult = index ($Nseq, $motif);			
		if ($indexresult > -1) {$ATG = $indexresult + $Nregion; print "found ATG at $ATG\n"}}
						
	if (($ATG < 0) && ($Repstart - 864 > 0)) {$ATG = $Repstart - 864} 
	elsif ($ATG < 0) { 
		if ($Repstart % 3 == 0) {$ATG = 3}
		elsif (($Repstart + 1) % 3 == 0) {$ATG = 2}
		elsif (($Repstart + 2) % 3 == 0) {$ATG = 1}
		};
	my $Nterm = substr $seq, $ATG, $Repstart - $ATG;
	
	
	#find C terminal
	
	my $STOP = 0;
	my @tempstops;
	my $Repend = $AoA[$j][$TALlen];
	my $Cseq =  substr $seq, $Repend, 1500;
	#print "$Nseq\n";

	foreach (@STOPmot) {
		my $motif = $_;
		my $indexresult = 0;
		my $offset = 0;
		while ($indexresult != -1 ) {
		$indexresult = index ($Cseq, $motif, $offset);		
		$offset = $indexresult + 1;
		if (($indexresult > -1) && ($indexresult % 3 == 0)) {push @tempstops, $indexresult + $Repend}
		}}


	if (@tempstops == 0) {$STOP = $Repend +  882} else {$STOP = min (@tempstops)}
		print "STOP found at $STOP\n";
	my $Cterm = substr $seq, $Repend, $STOP - $Repend;
	my $tot = $TALlen +0.5;

	my $spec; my $str;
	if (defined $species{$seqaccession}) {$spec = $species{$seqaccession}} else {$spec = $single_specie};
	if (defined $strains{$seqaccession}) {$str = $strains{$seqaccession}} else {$str = $single_strain};
	
	#my $nuheader = "$spec|$header"."_+_$ATG"."_$STOP"."_$tot";
	#my $nuheader = "$spec|$str|$header"."_"."$tot";
	#my $nuheader = "$header"."_+_$ATG"."_$STOP"."_$tot";
	my $nuheader = "$header";

	#print Nterm
	print $OUT1 ">$nuheader"." Nterminal\n$Nterm\n";

	#print Cterm
	print $OUT2 ">$nuheader"." Cterminal\n$Cterm\n";

	#print repeats nucleotide	

	print $OUT3 ">$nuheader"." Repeats\n"; foreach (@Repeat_array) {print $OUT3 "$_\n";}
	#print repeats aminoacid
	#MOD ALPQ, translate Nterm and Cterm
	my $AANTERM = translate_as_string ($Nterm);
	my $AACTERM = translate_as_string ($Cterm);
	my @NOLAST = @AArepeat;
	pop @NOLAST;

	print $OUT8 ">$nuheader"." Repeats\n"; print $OUT8 "$AANTERM\n"; foreach (@NOLAST) {print $OUT8 "$_\n";}; print $OUT8 "$AACTERM\n"; #MOD ALPQ NOV 2018
	
	my $whole = substr $seq, $ATG, $STOP-$ATG+3;
	print $OUT7 ">$nuheader\n$whole\n";
	my $wholeAA = translate_as_string ($whole);
	print $OUT6 ">$nuheader\n$wholeAA\n";
	#print RVDs
	
	print $OUT9 "$header\tRVDminer\tCDS\t$ATG\t$STOP\t.\t+\t0\tgene=TAL_$tot\n";

	print $OUT4 "$nuheader\t";

	for my $i ( 0 .. ($#RVD_array -1) ) {print $OUT4 "$RVD_array[$i]-";}
	print $OUT4 "$RVD_array[$#RVD_array]";
	if  (defined $species{$seqaccession}) {print $OUT4 "\t$species{$seqaccession}\t"} else {print $OUT4 "\t$single_specie\t"};
	if (defined $strains{$seqaccession}) {print $OUT4 "$strains{$seqaccession}\t"} else {print $OUT4 "$single_strain\t"};
	print $OUT4 "$seqaccession\t$seqgi\t@comments\n";

					}}

##########################################################################################################################################
for (my $j = 0; $j <= $#AoB; $j++) { #FOR EACH TAL in the MINUS strand
	my @comments; #comments to add about each sequence
	my $TALlen = $#{$AoB[$j]};
	if ($TALlen > 0) {
	print "TALlength = $TALlen\n";
	#extract repetitions and RVDs
	
	my @Repeat_len;
	my @Repeat_array;
	my @AArepeat;
	my @RVD_array;
	my $RVD;
	
	#find repeat length
	for ( my $i = 1 ; $i <= $TALlen ; $i++ ) {
		my $length = $AoB[$j][$i] - $AoB[$j][ $i - 1 ];
		push @Repeat_len, $length;
		
	}
	#find most common repetition length
		my(%count);
		foreach my $value (@Repeat_len) {
    		$count{$value}++;
		}
		my $common_len = (sort {$count{$b} <=> $count{$a}} @Repeat_len)[0];
		my $min_len = (sort {$a <=> $b} @Repeat_len)[0];
		my $max_len = (sort {$a <=> $b} @Repeat_len)[$#Repeat_len];
	

	for ( my $i = 0 ; $i <= $#Repeat_len ; $i++ ) {
		my $repeat = substr ($revcom, $AoB[$j][$i], $Repeat_len[$i]);
		push @Repeat_array, $repeat;
		my $aarepeat = translate_as_string ($repeat);
		push @AArepeat, $aarepeat;
		#print "$aarepeat\n";
		if ($Repeat_len[$i] == $common_len) {$RVD = substr ($aarepeat, 11, 2); push @RVD_array, $RVD;}
		#add asterisk to RVD if repeat is one aminoacid shorter
		elsif ($Repeat_len[$i] == $common_len - 3) {$RVD = substr ($aarepeat, 11, 1); push @RVD_array, $RVD."*"} 
		else {$RVD = substr ($aarepeat, 11, 2); push @RVD_array, $RVD; push @comments, "aberrant repeat;"; print "Aberrant or unrecognized repeat : $aarepeat\n"}
						}
	#get last repeat		
		my $last_repeat = substr ($revcom, $AoB[$j][$TALlen], $common_len/2);
		push @Repeat_array, $last_repeat;
		my $aalast_repeat = translate_as_string ($last_repeat);
		print "$aalast_repeat\n";
		push @AArepeat, $aalast_repeat;
		push @RVD_array, substr ($aalast_repeat, 11, 2);

	#find N terminal
	
	my $ATG = -1;
	my $Nregion;
	my $Nseq;

	my $Repstart = $AoB[$j][0];
	if ($Repstart - 1500 > 0) {$Nregion = $Repstart - 1500; $Nseq =  substr $revcom, $Nregion, 1500}
		else {$Nregion = 0; $Nseq = substr $revcom, 0, $Repstart}

	foreach (@Ntermmotifs) {
		my $motif = $_;
		my $indexresult = 0;
		$indexresult = index ($Nseq, $motif);			
		if ($indexresult > -1) {$ATG = $indexresult + $Nregion; print "found ATG at $ATG\n"}}
						
	
	if (($ATG < 0) && ($Repstart - 864 > 0)) {$ATG = $Repstart - 864} 
	elsif ($ATG < 0) { 
		if ($Repstart % 3 == 0) {$ATG = 3}
		elsif (($Repstart + 1) % 3 == 0) {$ATG = 2}
		elsif (($Repstart + 2) % 3 == 0) {$ATG = 1}
		};
	my $Nterm = substr $revcom, $ATG, $Repstart - $ATG;
		
	#find C terminal
	
	my $STOP = 0;
	my @tempstops;
	my $Repend = $AoB[$j][$TALlen];
	my $Cseq =  substr $revcom, $Repend, 1500;
	#print "$Nseq\n";

	foreach (@STOPmot) {
		my $motif = $_;
		my $indexresult = 0;
		my $offset = 0;
		while ($indexresult != -1 ) {
		$indexresult = index ($Cseq, $motif, $offset);		
		$offset = $indexresult + 1;
		if (($indexresult > -1) && ($indexresult % 3 == 0)) {push @tempstops, $indexresult + $Repend}
		}}

	if (@tempstops == 0) {$STOP = $Repend +  882} else {$STOP = min (@tempstops)}
		print "STOP found at $STOP\n";
		my $Cterm = substr $revcom, $Repend, $STOP - $Repend;
	my $tot = $TALlen +0.5;
	
	my $minATG = $seqlen - $ATG;
	my $minSTOP = $seqlen - $STOP;

	
	my $spec; my $str;
	if (defined $species{$seqaccession}) {$spec = $species{$seqaccession}} else {$spec = $single_specie};
	if (defined $strains{$seqaccession}) {$str = $strains{$seqaccession}} else {$str = $single_strain};
	
	#my $nuheader = "$spec|$str|$header"."_"."$tot";
	
	#my $nuheader = "$spec|$header"."_-_$minSTOP"."_$minATG"."_$tot";
	#my $nuheader = "$header"."_-_$minSTOP"."_$minATG"."_$tot";
	my $nuheader = "$header";

	#print Nterm
	print $OUT1 ">$nuheader"." Nterminal\n$Nterm\n";

	#print Cterm
	print $OUT2 ">$nuheader"." Cterminal\n$Cterm\n";

	#print repeats nucleotide	

	print $OUT3 ">$nuheader"." Repeats\n"; foreach (@Repeat_array) {print $OUT3 "$_\n";}
	#print repeats aminoacid
	my $AANTERM = translate_as_string ($Nterm);
	my $AACTERM = translate_as_string ($Cterm);
	my @NOLAST = @AArepeat;
	pop @NOLAST;

	print $OUT8 ">$nuheader"." Repeats\n"; print $OUT8 "$AANTERM\n"; foreach (@NOLAST) {print $OUT8 "$_\n";} ; print $OUT8 "$AACTERM\n"; #MOD ALPQ NOV 2018
	#print $OUT8 ">$nuheader"." Repeats\n"; foreach (@AArepeat) {print $OUT8 "$_\n";}	
	
	my $whole = substr $revcom, $ATG, $STOP-$ATG+3;
	print $OUT7 ">$nuheader\n$whole\n";
	my $wholeAA = translate_as_string ($whole);
	print $OUT6 ">$nuheader\n$wholeAA\n";
	#print RVDs
	
	print $OUT4 "$nuheader\t";

	print $OUT9 "$header\tRVDminer\tCDS\t$minSTOP\t$minATG\t.\t-\t0\tgene=TAL_$tot\n";

	for my $i ( 0 .. ($#RVD_array -1) ) {print $OUT4 "$RVD_array[$i]-";}
	print $OUT4 "$RVD_array[$#RVD_array]";
	if  (defined $species{$seqaccession}) {print $OUT4 "\t$species{$seqaccession}\t"} else {print $OUT4 "\t$single_specie\t"};
	if (defined $strains{$seqaccession}) {print $OUT4 "$strains{$seqaccession}\t"} else {print $OUT4 "$single_strain\t"};
	print $OUT4 "$seqaccession\t$seqgi\t@comments\n";

					}}

}
	#print "$aaseq\n";
	#if ($header =~ /^gi/) {	#for genbank protein fastas
	#$seqgi = (split(/\|/,$id))[1];
	#$seqaccession = (split (/\./, (split(/\|/,$id))[3]))[0];
	#if (not defined $species{$seqaccession}) {$species{$seqaccession} = (split(/[\[\]]/,$header))[1]};
	#			}

my $cmd = "seqret -sequence $ARGV[0] -feature -fformat gff -fopenfile $ARGV[1]_TALs.gff3 -osformat genbank -outseq $ARGV[1].gbk";
system ("$cmd");


	
