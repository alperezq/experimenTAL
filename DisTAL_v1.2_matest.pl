#!/usr/bin/perl -w
use warnings;
use strict;
use Getopt::Std;
use Statistics::R;
use List::MoreUtils qw(uniq);
use List::Util qw( min max );
use Algorithm::NeedlemanWunsch;
use Bio::Perl;
use POSIX qw(ceil);

##################################################################
# DisTAL							 #
# by Alvaro L Perez-Quintero                       		 #
#   < alperezqui at gmail dot com >                		 #
#   Created							 #
#   Last modification August 29 2016	     			 #
#								 #
#This script uses the program  ARLEM  version 1.0    		 #
# by Mohamed I. Abouelhoda  (C) 2007      			 #
# available at http://www.nubios.nileu.edu.eg/tools/wami/	 #
##################################################################


sub printCMD() {
print STDERR "\n

Usage: DisTAL [options] TALfile Outname(default=output) AdditionalInfo(optional)

############################################################################
			QueTAL - DisTALv1.2: 

This script takes a fasta file of either nucleotide or aminoacid sequences of
TAL effectors and calculates distances based on the repeats of these TALs
by considering each repeat as a unit. This version also considers the N and C temrinal region
as repeats if they are available

This script uses the program  ARLEM  version 1.0    				
by Mohamed I. Abouelhoda  (C) 2007      			 		
available at http://www.nubios.nileu.edu.eg/tools/wami/

Requires: Perl modules Statistics::R, Algorithm::NeedlemanWunsch; Bio::Perl;
	 R, R library 'ape'
#############################################################################

TALfile: Fasta file with either nucleotide or aminoacid sequences of TAL effectors
AdditionalInfo: Tab separated file with fasta Ids in the first column and information 
to be used to color the outpu tree (e.g. specie or strain)

-c :[T or F]: compare the input TALs against a set of available TAL sequences from the ncbi database (default = F)
-l :[number]: if d=T, number of target TAls form database to output (default = 5)
-n :[f, d, u, p or c] Tree format: organze the output tree in a fan (f), phylogra, (p), dendrogram (d), unrooted (u)(default) or cladogram (c) format
-r :[T or F]: use RVDs; RVDs could be excluded from the analyses to avoid the selection effects on these aminoacids (default = T)


Outputs 

Outname.matrix: Distance matrix based on TAL repeats
Outname.tre: nj tree in newick format
Outname.pdf: nj tree in pdf format
Outname.hits (if -d = T): Similar TALs found in the database according to specificities 
Outname_Repeatscode: Unique repeats found in input file asigned to a uniq number
Outname_CodedRepeats: Original sequence coded as series of unique repeats according to repatcode

";
exit;
}

if (@ARGV < 1) {
	printCMD();
}


#######Get arguments

my %options=();
getopts("hc:n:m:l:i:d:r:", \%options);

my $newMatrix = '';
my $database = '';
my $Treetype = '';
my $lim = '';
my $indel = '';
my $dup = '';
my $ruse = '';
my $indefault = 10;
my $dupdefault = 10;


if (defined $options{c}) {$database = $options{c}} else {$database = "F"};
if (defined $options{n}) {$Treetype  = $options{n}} else {$Treetype = "u"};
if (defined $options{m}) {$newMatrix = $options{m}} else {$newMatrix = "F"};
if (defined $options{l}) {$lim = $options{l}} else {$lim = 10};
if (defined $options{i}) {$indel = $options{i}} else {$indel = $indefault};
if (defined $options{d}) {$dup = $options{d}} else {$dup = $dupdefault};
if (defined $options{r}) {$ruse = $options{r}} else {$ruse = "T"};

my $Outdir = "Outputs";
my $Outname;
if ( -d "$Outdir") {} else { system ("mkdir $Outdir")};
if ($ARGV[1]) {$Outname = "$Outdir/$ARGV[1]"} else {$Outname = "$Outdir/Output"};

##############################Read input file and format as needed

my @RAW;
open my $IN0, "<$ARGV[0]" or die "Can not open $ARGV[0], file doesn't exist\n";

while ( my $line = <$IN0> ) 	{
	chomp $line;
	push @RAW, $line;
		}




###Detect input format

my $countaa = 0;
my $countnuc = 0;
my $countfas = 0; 
for (my $i = 0; $i <= 10; $i++) {
 if ((defined $RAW[$i])  and ($RAW[$i] =~ '>')) {$countfas++;}
	elsif ((defined $RAW[$i]) and ("$RAW[$i]" =~ /[HDLRP]+/)) {$countaa++}
	elsif ((defined $RAW[$i]) and ("$RAW[$i]" =~ /[ATCG]+/i)) {$countnuc++}

}


if (($countfas == 0) or ($countnuc == 0 && $countaa == 0)) {die "Input file is neither an aminoacid or nucleotide fasta\n";}

#convert to aminoacid if necessary

my $AAref; my $IDsref;

if ($countaa <= $countnuc) {
print STDERR "Input file is a nucleotide fasta, converting to aminoacid\n";
($AAref, $IDsref) = read_fasta(\@RAW, "translate");
}#if nucleotide

else 

{
($AAref, $IDsref) = read_fasta(\@RAW);
}

my %AA = %$AAref;
my @rIDs = @$IDsref;

my @Allrepeats;
my @lines;

open my $Inital, "<Info/Initial_TAL_strings.txt";
my @Ini_repeat;

while ( my $line = <$Inital> ) {
	chomp $line;
	push @Ini_repeat, $line;}



####################FORMAT REPEATS

for ( my $i = 0 ; $i <= $#rIDs ; $i++ ) {
my $id = "$rIDs[$i]";
my $seq = $AA{$rIDs[$i]};

my @position_array; #array that will store all the postitions where a string for a TAL repetition is found
#look for each possible TAL repetition initiation
foreach (@Ini_repeat) {
		my $ini_motif = $_;
		my $indexresult = 0;
		my $offset = 0;
  #use index function in a loop to find all occurences of the TAL initial string
		while ( $indexresult != -1 ) {
			$indexresult = index( $seq, $ini_motif, $offset );			
			$offset = $indexresult + 1;			
			if ( $indexresult > -1 ) { push @position_array, $indexresult }
		}
	}
  #create array with uniq positions where the the initial TAL string was found
	my @POS = uniq( sort { $a <=> $b } (@position_array) );

  #get repeats lengths
	my @Repeat_len;

	for ( my $i = 1 ; $i <= $#POS ; $i++ ) {
		my $length = $POS[$i] - $POS[ $i - 1 ];
		push @Repeat_len, $length;
		
	}
  #find most common repetition length
		my(%count);
		foreach my $value (@Repeat_len) {
    		$count{$value}++;
		}
		my $common_len = (sort {$count{$b} <=> $count{$a}} @Repeat_len)[0];
		my $min_len = min @Repeat_len;
		my $max_len = max @Repeat_len;

   #exclude sequences for which no TAL repetition was found 
		if (not defined $min_len) 
	{print "no TAL repetitions were found for sequence $id, it was excluded from further analyses\n"}
		else {

	my @Repeat_array;
	#Add Nterminal to repeat array
	
	if ("$POS[0]" > 0){
	my $Nter = substr ($seq, 0, $POS[0]);
	print "$Nter\n";
	push @Repeat_array, $Nter;
	}

	for ( my $i = 0 ; $i <= $#Repeat_len ; $i++ ) {
		my $times = ($Repeat_len[$i] / $common_len);
		my $repeat;

		if ($times > 1.5) { #deal with XXXX repeats
		my $ini = "$POS[$i]";
		for  (my $k = 1 ; $k <= $times ; $k++ ) {
						if ($ruse eq "T") {$repeat = substr ($seq, $ini, $common_len)} 
						elsif ($ruse eq "F") {$repeat = substr ($seq, $ini, ($common_len - 2))}; 
						$ini = $ini + ($common_len -2);
						push @Repeat_array, $repeat;
							}}
		else { #no XXX repeats
		
			if ($Repeat_len[$i] == $common_len - 1) { #repats with N* or H*
				if ($ruse eq "T") {$repeat = substr ($seq, $POS[$i], $Repeat_len[$i])}
				elsif ($ruse eq "F") { #avoid RVDs
				my $rep1 = substr ($seq, $POS[$i], 11);
				my $rep2 = substr ($seq, $POS[$i] + 12, ($Repeat_len[$i]-12));
				$repeat = join ("", $rep1, $rep2);
							}
				push @Repeat_array, $repeat; #print "$repeat\n"
				}

			else {
				if ($ruse eq "T") {$repeat = substr ($seq, $POS[$i], $Repeat_len[$i])}
				elsif ($ruse eq "F") { #avoid RVDs		
				my $rep1 = substr ($seq, $POS[$i], 11);
				my $rep2 = substr ($seq, $POS[$i] + 13, ($Repeat_len[$i]-13));
				$repeat = join ("", $rep1, $rep2);
							}
				push @Repeat_array, $repeat; #print "$repeat\n" 
				}

				}
				}
	#get last repeat		
		my $last_repeat;
		if ($ruse eq "T") {$last_repeat = substr ($seq, $POS[$#POS], $common_len)}
		elsif ($ruse eq "F") { #avoid RVDs		
				my $rep1 = substr ($seq, $POS[$#POS], 11);
				my $rep2 = substr ($seq, $POS[$#POS] + 13, ($common_len-11));
				$last_repeat = join ("", $rep1, $rep2);
							}
		push @Repeat_array, $last_repeat;

###ADD Cterminal

	my $Cter = substr ($seq, $POS[$#POS]+$common_len);
	print "$Cter\n";
	push @Repeat_array, $Cter;


	#print repeats
	if ( ($max_len == $common_len) and ($min_len > $common_len -2)) 
		{push @lines, "$id"; foreach (@Repeat_array) {push @lines, "$_"; push @Allrepeats, "$_"}} 
	else 	{#print "Longer or unusual repeats were found for sequence $id, it was still kept for analyses\n";
		push @lines, "$id"; foreach (@Repeat_array) {push @lines, "$_"; push @Allrepeats, "$_"}}

	#print "$id\n"; foreach (@Repeat_array) {print "$_\n"};
				}



}

#Read additional Info file

my $Addinfoflag = "F";
my @Addinfo;

if (defined $ARGV[2] && ($database eq "F")) {
my $count = 0;
open my $ADD, "<$ARGV[2]" or die "Can not open $ARGV[0], file doesn't exist\n";
while ( my $line = <$ADD> ) {
	chomp $line;
	my $addID = (split (/\t/, $line))[0];
	if ($addID  =~ '>') {} else {$addID = ">$addID"}
	if (exists $AA{$addID}) {
	push @Addinfo, (split (/\t/, $line))[1];
	$count ++;
		}

}
if ($count == keys %AA) {$Addinfoflag = "T";} else {print STDERR "\tIDs in additional info file do no match input fasta file, additional info was no used\n"}
}


#ADD PUBLIC TALS IF DATABASE OPTION IS TRUE

my $HITS;

if ($database ne "F") {
open $HITS, ">$Outname.hits" or die "Can not open $Outname.hits, file doesn't exist\n";
print STDERR "\t Including TALs from public databases in the analyses...\n";
open my $PUBLIC, "<Info/Public_TALs_repeats_AA.FASTA";
while ( my $line = <$PUBLIC> ) {
	chomp $line;
	push @lines, $line;
	if ($line !~ '>'){	
	push @Allrepeats, $line;

		}
}
}

###############Get uniqRVDs


my @Uniqrepeats = uniq( sort (@Allrepeats) ); #get uniq RVDs

my $MATfile;
my $CLASSINfile;
my $CLASSfileUsed = "$Outname"."_Repeatscode.txt";  
my %ClassRep;

#####If a new distance matrix for input repeats is required create it

if ($newMatrix ne "F") {

print STDERR "\t Creating distance matrix for repeats...\n";

$MATfile = "$Outname"."_Repeatmatrix.mat";
open my $OUT0, ">$MATfile" or die "Can not open $MATfile, file doesn't exist\n";
open my $OUT1, ">$CLASSfileUsed" or die "Can not open $CLASSfileUsed, file doesn't exist\n";

my $num = $#Uniqrepeats + 1;
print $OUT0 "# Type no. $num\n# Types ";
for (my $i = 0; $i <= $#Uniqrepeats; $i++) {print $OUT0 "$i "; 
					print $OUT1 "$i\t$Uniqrepeats[$i]\n"; #print used seqs
					$ClassRep{$Uniqrepeats[$i]} = $i;
					};
print $OUT0 "\n# Indel align $indel
# Indel hist $indel
# Dup $dup
# matrix 
";

for (my $i = 0; $i <= $#Uniqrepeats -1; $i++)		{
my $seqa = "$Uniqrepeats[$i]";
for (my $j = $i + 1; $j <= $#Uniqrepeats; $j++)		{
my $seqb = "$Uniqrepeats[$j]";
#print "$i._.$j\t";
my @Sa = split //, $seqa; # make an array
my @Sb = split //, $seqb;

my $matcher = Algorithm::NeedlemanWunsch->new(\&score_sub);
my $score = $matcher-> align(
                \@Sa,
                \@Sb,
		);

my $maxlen = (sort { $a <=> $b } ($#Sa, $#Sb))[-1] + 1;
#my $normscore = int((($maxlen - $score) / $maxlen) * 100 )+0.5); #mod
my $normscore = (($maxlen - $score) / $maxlen) * 100 ; #mod
#$normscore = sprintf("%.0f", $normscore);
$normscore = ceil($normscore);
print $OUT0 "$normscore "; 

}
print $OUT0 "\n";
}
close $OUT0;
}

###If no new matrix is required match the input repeats to existing one from the databases

else {

print STDERR "\t Comparing repeats to existing distance matrix...\n";

#$MATfile = "Info/ALLMatrixaa";
#$CLASSINfile = "Info/ALLUniqaa";

#$MATfile = "Info/XOMATRIX.mat"; #MOD ALPQ FOR NU PACBIO!!!!!
#$CLASSINfile = "Info/XOSEQS.txt";

$MATfile = "Info/XOLS_Repeatmatrix.mat"; #MOD AGAIN FOR AVRXA7, CHANGE BACK!!
$CLASSINfile = "Info/XOLS_Repeatscode.txt";


open my $IN1, "<$CLASSINfile"; #Read Uniq repeats from database

open my $OUT2, ">$CLASSfileUsed";

my @Repseq;

while (my $line = <$IN1>){
	chomp $line;
	my @tmp = split (/\t/, $line);
	$ClassRep{$tmp[1]} = $tmp[0]; #save uniq sequences and their ID in a hash
	push @Repseq, $tmp[1];
				}
close $IN1;

my @NewRep;

my $max = $#Repseq;

foreach (@Uniqrepeats) { #check if there are new repeats
my $Rep = "$_";
if (defined $ClassRep{$Rep}) {print $OUT2 "$ClassRep{$Rep}\t$Rep\n"}
			else {
			push @NewRep, $Rep;
			$max++;
			$ClassRep{$Rep} = $max;
			print $OUT2 "$ClassRep{$Rep}\t$Rep\n"
					}
}
#####if new repeats were found, create new distance matrix by aligning the new repeats against the existing repeats

if ($#NewRep >= 0) {

print STDERR "\t New repeats found, calculating distances for repeats and creating new matrix...\n";


open my $IN3, "<$MATfile";
$MATfile = "$Outname"."_Newmatrix.mat";
open my $OUT3, ">$MATfile";

my @OryMat;
while ( my $line = <$IN3> ) {chomp $line; if ($line !~ '#') {push @OryMat, $line};}

my $num = $#Repseq + $#NewRep + 2;
print $OUT3 "# Type no. $num\n# Types ";
for (my $i = 0; $i <= $num - 1; $i++) {print $OUT3 "$i ";};
print $OUT3 "\n# Indel align $indel
# Indel hist $indel
# Dup $dup
# matrix 
";


for (my $z = 0; $z <= $#NewRep; $z++) {

print STDERR "Adding Rep Number $z out of $#NewRep\n";
my $seqb = $NewRep[$z];

for (my $i = 0; $i <= $#Repseq; $i++) {
	my $seqa = $Repseq[$i];  #########MAKE ALIGMENT BETWEEN SEQA AND SEQB 

my @Sa = split //, $seqa; # make an array
my @Sb = split //, $seqb;

my $matcher = Algorithm::NeedlemanWunsch->new(\&score_sub);
my $score = $matcher-> align(
                \@Sa,
                \@Sb,
		);

my $maxlen = (sort { $a <=> $b } ($#Sa, $#Sb))[-1] + 1;
my $normscore = (($maxlen - $score) / $maxlen) * 100 ; #mod
$normscore = ceil($normscore);

if (defined $OryMat[$i]) {$OryMat[$i] = "$OryMat[$i]$normscore "} else {$OryMat[$i] = "$normscore "}
		}
push @Repseq, $seqb;
}
foreach (@OryMat) {print $OUT3 "$_\n";}
}#new reps found

elsif ($indel != $indefault or $dup != $dupdefault)
{

open my $IN3, "<$MATfile";
$MATfile = "$Outname"."_Newmatrix.mat";
open my $OUT3, ">$MATfile";
my @OryMat;
while ( my $line = <$IN3> ) {chomp $line; if ($line !~ '#') {push @OryMat, $line};}
my $num = $#Repseq + $#NewRep + 2;
print $OUT3 "# Type no. $num\n# Types ";
for (my $i = 0; $i <= $num - 1; $i++) {print $OUT3 "$i ";};
print $OUT3 "\n# Indel align $indel
# Indel hist $indel
# Dup $dup
# matrix 
";
foreach (@OryMat){print $OUT3 "$_\n"}
}


}#else no newmatrix

############# take original sequences and transform them into coded alphabet


#for (my $i = 0; $i <= $#Uniqrepeats; $i++) {$ClassRep{$Uniqrepeats[$i]} = $i}; MOD ALPQ !!!! ERRor when comparing to existing Databases
print STDERR "\t Coding TAL repeats into custom alphabet...\n";

my $CodedRepeats = "$Outname"."_CodedRepeats.fa";
open my $OUT4, ">$CodedRepeats";
my @IDs;
my %TALlen;
my $idcount = -1;

my $zero = "$lines[0]";
print $OUT4 "$zero\n"; $zero =~ s/>//g; push @IDs, $zero; $idcount++; #first line

for (my $i=1;$i<=$#lines;$i++) {
my $line = $lines[$i];
if ($line =~ '>'){
	print $OUT4 "\n$line\n"; $line =~ s/>//g; push @IDs, $line; $idcount++;
		}
else {print $OUT4 "$ClassRep{$line} "; $TALlen{$IDs[$idcount]}++} #get repeat length
}

#for (my $i=0;$i<@IDs;$i++) {print  "$IDs[$i]\t$TALlen{$IDs[$i]}\n"}



######Run Arlem and process output
print STDERR "\t Running alignment of coded repeats using Arlem v.1 by Mohamed I. Abouelhoda (C) 2007 ...\n";

my @Arlemoutput = `./arlem -f $CodedRepeats -cfile $MATfile -align -insert -showalign`;

my $Alignmatrix = "$Outname.mat";
open my $OUT5, ">$Alignmatrix";

my @Scores;
my @lookformax;

foreach (@Arlemoutput){
my $res = $_;
print "$res\n";
if ($res =~ "Score of aligning") {
my $s1 = (split(/:/,$res))[1];
my $s2 = (split(/:/,$res))[2];
my $N1 = (split(/,/,$s1))[0];
my $N2 = (split(/ /,$s2))[0];
$Scores[$N1][$N2] = (split(/=/,$res))[1];
$Scores[$N1][$N2] =~ s/\n//g;
push @lookformax, "$Scores[$N1][$N2]";
}}

my $max = max (@lookformax);

#############create aligment matrix for trees

print STDERR "\t Creating TAL distance matrix ...\n";
print $OUT5 "matrix\t";

for (my $i=0;$i<@IDs;$i++) {
	print $OUT5 "$IDs[$i]\t";
}
	print $OUT5 "\n";
for (my $i=0;$i<@IDs;$i++) {
	my $n = ">$IDs[$i]";
	my %scores;
	print $OUT5 "$IDs[$i]\t";
	for (my $j=0;$j<@IDs;$j++) {
	my $o = ">$IDs[$j]";
	my $val;
	if ($i == $j) {print $OUT5 "0\t"; $val = "0"}
	my $maxTALlen = max ($TALlen{$IDs[$i]}, $TALlen{$IDs[$j]}); #look for TALlength of paired TALs to normalize score
	#print "$maxTALlen\n";
	#if (defined $Scores[$i][$j]) {$val = $Scores[$i][$j] / $max; print $OUT5 "$val\t"} #nromalize by maximum distance in matrix
	#elsif (defined $Scores[$j][$i]) {$val = $Scores[$j][$i] / $max; print $OUT5 "$val\t"} #modified

	if (defined $Scores[$i][$j]) {$val = $Scores[$i][$j] / $maxTALlen ; print $OUT5 "$val\t"} #normalized by length
	elsif (defined $Scores[$j][$i]) {$val = $Scores[$j][$i]/ $maxTALlen ; print $OUT5 "$val\t"}

	#if (defined $Scores[$i][$j]) {$val = $Scores[$i][$j]; print $OUT5 "$val\t"} #not normalized
	#elsif (defined $Scores[$j][$i]) {$val = $Scores[$j][$i]; print $OUT5 "$val\t"}


	if ($database ne "F" && exists $AA{$n} && not exists $AA{$o}) {$scores{$o} = $val;}
				}

	if ($database ne "F" && exists $AA{$n}) { #find closest TALs in database
			print $HITS "The TALs more similar repeats to $IDs[$i] in our database are:\nID\tDistance\n";
			my @SimilarID; my @SimilarScore;
			foreach my $name (sort { $scores{$a} <=> $scores{$b} } keys %scores) {push @SimilarID, "$name"; push @SimilarScore, "$scores{$name}";}
			for (my $j=0;$j<=$lim - 1;$j++)
						{
				print $HITS "$SimilarID[$j]\t$SimilarScore[$j]\n";
						}
						}
	print $OUT5 "\n";
				}



#####BUILD TREE using R

print STDERR "\t Bulding tree...\n";

my $R = Statistics::R->new(); #bridge to use with R




my $cmds1 = <<EOF;
library("ape")
exprs <- as.matrix(read.table("$Alignmatrix", header=TRUE, sep = "\t",row.names = 1, as.is=TRUE))
tr <- nj((exprs))
#bt<-boot.phylo(tr, exprs, FUN = function(xx) nj(exprs), B = 1000)
write.tree(tr, file="$Outname.tre")
EOF

my $run1 = $R->run($cmds1);

if ($Addinfoflag eq "T") { #Use additional info to color the tree

$R->set('labels', \@Addinfo);

my $cmds2 = <<EOF;

flabels<-as.factor(labels)
colvector<-rainbow (nlevels(flabels))[as.integer(flabels)]
pdf("$Outname.pdf")
plot(tr,"$Treetype", tip.color=colvector, show.tip.label=TRUE, cex=0.4)
tiplabels(pch = 21, bg=colvector, col = "black", lwd = 0.3, cex=1.2)
add.scale.bar()
#nodelabels(bt,frame = "none", cex = 0.5, col = "blue")
legend(x = 'bottomright',legend = as.character(levels(flabels)),col = rainbow (nlevels(flabels)), pch = par("pch"), bty = 'n', xjust = 1, cex = 0.7)
dev.off()
EOF

my $run2 = $R->run($cmds2);
}


else {

my $cmds2 = <<EOF;
pdf("$Outname.pdf")
plot(tr,"$Treetype", show.tip.label=TRUE, cex=1, use.edge.length = TRUE)
add.scale.bar()
#nodelabels(bt,frame = "none", cex = 0.5, col = "blue")
dev.off()
EOF
my $run2 = $R->run($cmds2);
}




### scoring function:
sub score_sub {
  if (!@_) {
    return 0;            # gap penalty
  }
  ## mismatch scores -1, match +1
  return ($_[0] eq $_[1]) ? 1 : -1; #dont change!
} 


###reading fasta



sub read_fasta { #takes two arguments; one fasta array and a flag to translate or not to protein
my $inlines = shift;
my @lin= @{$inlines};
my @rIDs;
my $trans = shift;
my %fast;
for (my $i = 0; $i <= $#lin - 1; $i++)
			{
if ($lin[$i] =~ '>') {#process fasta
	push @rIDs, "$lin[$i]";
my $j = $i + 1;
my $string = '';
do
				{
		if ((defined $lin[$j]) && ($lin[$j] !~ '>')) {
			$string = join ("", $string, "$lin[$j]");
				$j++;
					}
else {$j= 0}
		}while($j != 0);
if (defined $trans) {my $prot_obj; $prot_obj = translate_as_string ($string); $fast{$lin[$i]} = "$prot_obj";}
else {$fast{$lin[$i]} = "$string"};
}
}
return (\%fast, \@rIDs);
}
