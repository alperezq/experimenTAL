#!/usr/bin/perl -w
use warnings;
use strict;
use List::UtilsBy qw(max_by);
use List::MoreUtils qw(first_index);

#usage perl Analyze_TALs2.pl Aligned_TEV_ALL Coded_Reps_withgroups_4.8.fa";

#Store groups 

open my $IN0, "<TALgroups.txt";

my @RAW;


while ( my $line = <$IN0> ) {
	chomp $line;
	push @RAW, $line;
}
close $IN0;

my %TAL_Group; #group for each tal ID
my %Group_TAL; #arrays of Ids for each group

#Store hash of groups

for (my $i =1; $i <= $#RAW; $i++) {#mod ALPQ
	my $line = $RAW[$i];
	my @tmp = split (/\s/, $line);
	my $group = "G$tmp[$#tmp]";
	my $nexttolast = $#tmp -1;
	#my $reid = join (" ",@tmp[0..$nexttolast]);	;#MOD ALPQ, error for Nupacbio with some strains having spaces in their names
	my $reid = "$tmp[0]";
	$TAL_Group{$reid} = $group;
	push @{ $Group_TAL{$group} }, "$reid|$group";
	#$Group_counts{$group}++;
}


#Store hash of alignements

open my $IN1, "<$ARGV[0]";

my %ALIGNS;

while ( my $line = <$IN1> ) {
	chomp $line;
	my @tmp = split (/\t+/, $line);
	$ALIGNS{$tmp[0]}{$tmp[1]} = "$tmp[2]\t$tmp[3]";
	$ALIGNS{$tmp[1]}{$tmp[0]} = "$tmp[3]\t$tmp[2]";
}

close $IN1;

#Store hash of unaligned sequences


my %SEQS;

open my $IN2, "<$ARGV[1]";
my @FRAW;

while ( my $line = <$IN2> ) {
	chomp $line;
	push @FRAW, $line;
}

close $IN2;

for (my $i =0; $i <= $#FRAW; $i++) {#mod ALPQ
	my $line = $FRAW[$i];
	my $ID;
	if ($line =~ ">") {$ID = $line; $ID =~ s/\>//g, $SEQS{$ID} = "$FRAW[$i+1]";}
}

open my $OUT1, ">Multiple_align";

#for each group with more than one member find the longest alignement

foreach my $Group (sort keys %Group_TAL) {
		my @Garr =  @{ $Group_TAL{$Group} };
	if ($#Garr == 0) {print $OUT1 "Group = $Group\n$Garr[0]   $SEQS{$Garr[0]}\n"} 
	else {

			#find longest alignements for each group
			
			my %LENs;
			for (my $i =0; $i <= $#Garr ; $i++) {
			my $seqa= $Garr[$i];
			#$LENs{$seqa} = 0;
			for (my $j = 0; $j <= $#Garr; $j++) {
			if ($i != $j) {
			my $seqb= $Garr[$j];
			my @tmpa = split (/\t/, $ALIGNS{$seqa}{$seqb});
			my @tmpb = split (/\s+/, $tmpa[1]);
			my @tmpc = split (/\s+/, $tmpa[0]);
			my $lenb = $#tmpb;
			my $lena = $#tmpc;
			$LENs{"$seqa"."&$seqb"} = $lenb + $lena;

			}

			}}

			my $maxid = max_by { $LENs{$_} } keys %LENs;
			my $max1 = (split (/&/, $maxid))[0];
			my $max2 = (split (/&/, $maxid))[1];
			

			#my $maxind = say first_index { $_ eq "$maxid" } @Garr;
			
			#find longest ID
			my $len = length "$Garr[0]";
				for my $str (@Garr) {
    			if ( length($str) > $len ) {
        		$len = length($str);
    			}
			}
			my $l= $len;

			print $OUT1 "Group = $Group\n";
			my @IN;
			for (my $i =0; $i <= $#Garr; $i++) {if (($Garr[$i] ne $max1) && ($Garr[$i] ne $max2)){push @IN, $i}}

			print "Maximum for group $Group are $max1 and $max2, total members are $#Garr, non max members are $#IN\n ";
			
			#choose longest between 2
			my @seql1 = split (/\s+/, $SEQS{$max1});
			my @seql2 = split (/\s+/, $SEQS{$max2});
			
			my $top; my $bot;

			if ($#seql1 > $#seql2) {$top = $max1; $bot = $max2} else {$top = $max2; $bot = $max1};

			#Print alignements
			my $al0 = $ALIGNS{$top}{$bot}; #print longest alignment
			my @tmpal0 = split (/\t/, $al0);

			printf $OUT1 "%-${l}s",  $top;
			my @arr1 = split (/\s/, $tmpal0[0]);
			printf $OUT1 "%5s", $_  for @arr1;
			print $OUT1 "\n";

			printf $OUT1 "%-${l}s", $bot;
			my @arr2 = split (/\s/, $tmpal0[1]);
			printf $OUT1 "%5s", $_  for @arr2;
			print $OUT1 "\n";
			
			if ($#IN >= 0) {
			for (my $i =0; $i <= $#IN; $i++) {
			my $p = $IN[$i];

			my $al1 = $ALIGNS{$top}{$Garr[$p]}; #print other alignements, find longest alignement to top or bottom
			my $al2 = $ALIGNS{$bot}{$Garr[$p]};
			my @tmpal1 = split (/\t/, $al1);
			my @tmpal2 = split (/\t/, $al2);
			my @topc = 	split (/\s/, $tmpal1[1]); 
			my @botc =  split (/\s/, $tmpal2[1]);
			my @arr;
			if ($#topc >= $#botc) {@arr = @topc} else {@arr = @botc};
			printf $OUT1 "%-${l}s", $Garr[$p];
			#my @arr = split (/\s/, $tmpal[1]);
			printf $OUT1 "%5s", $_  for @arr;
			print $OUT1 "\n";

			}#if group has more than 2
			}#for each group member
	}#else group with more than 1
}

