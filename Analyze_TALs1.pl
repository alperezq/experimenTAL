#!/usr/bin/perl -w
use warnings;
use strict;
#use Bio::Perl;
#use Bio::SeqIO;
#use Bio::Seq;
#use List::MoreUtils qw(uniq);
#use Algorithm::NeedlemanWunsch;
#use POSIX qw(ceil);

use Statistics::R;

#usage perl Analyze_TALs1.pl TEV_ALL_1.mat 5 TEV_ALL_1_CodedRepeats.fa F


## assign groups to files

my $MAT = "$ARGV[0]";
my $Nam = (split(/\./,"$MAT"))[0];
my $height = "$ARGV[1]";

my $R = Statistics::R->new(); #bridge to use with R

my $cmds1 = <<EOF;

library(reshape2)

TEV_ALL_1 <- as.matrix(read.table("$MAT", header=TRUE, sep = "\t",row.names = 1, as.is=TRUE))
TEV_ALL_1 <- TEV_ALL_1[1:ncol(TEV_ALL_1)-1,]
a <-(hclust(as.dist(TEV_ALL_1)))
CUT<-as.data.frame(cutree(a,h="$height"))
Groups <- as.data.frame(cbind(rownames(TEV_ALL_1),CUT[,1]))
colnames(Groups)<-c("TAL","Group")
write.table(Groups, "TALgroups.txt", row.names = FALSE,quote = FALSE)
EOF

my $run1 = $R->run($cmds1);

#Store groups 

open my $IN0, "<TALgroups.txt";

my @RAW;


while ( my $line = <$IN0> ) {
	chomp $line;
	push @RAW, $line;
}

my %TAL_Group;

#Store hash of groups

for (my $i =1; $i <= $#RAW; $i++) {#mod ALPQ
	my $line = $RAW[$i];
	my @tmp = split (/\s/, $line);
	$TAL_Group{"$tmp[0]"} = "G$tmp[1]";
}

#Assign groups to coded repeats 


open my $IN1, "<$ARGV[2]";
my @FRAW;

while ( my $line = <$IN1> ) {
	chomp $line;
	push @FRAW, $line;
}


open my $OUT0, ">Coded_Reps_withgroups_$height.fa";
open my $OUT1, ">Coded_Reps_withgroups_nolast_$height.fa";


#print files with group codes

for (my $i =0; $i <= $#FRAW; $i++) {#mod ALPQ
	my $line = $FRAW[$i];
	my $ID;
	if ($line =~ ">") {$ID = $line; $ID =~ s/\>//g; print $OUT0 ">$ID|$TAL_Group{$ID}\n"; print $OUT1 ">$ID|$TAL_Group{$ID}\n"}
	else {
		print $OUT0 "$line\n";
		my @tmp = split (/\s/, $line);
		for (my $j =0; $j <= $#tmp -1; $j++) 
			{print $OUT1 "$tmp[$j] "}
		print $OUT1 "\n";
		}
}

my $cfile = "$Nam"."_Repeatmatrix.mat";

print "Still going\n";

#Run arlem and align
if ("$ARGV[3]" eq "T") {
system ("./arlem -f Coded_Reps_withgroups_$height.fa -cfile $cfile -align -insert -showalign >TEV_ALL_arlout");
system ("perl pp_map_result_ALPQ.pl TEV_ALL_arlout Coded_Reps_withgroups_$height.fa >Aligned_TEV_ALL");
#system ("./arlem -f Coded_Reps_withgroups_nolast_$height.fa -cfile $cfile -align -insert -showalign >TEV_ALL_arlout");
#system ("perl pp_map_result_ALPQ.pl TEV_ALL_arlout Coded_Reps_withgroups_nolast_$height.fa >Aligned_TEV_ALL_nolasts");
}

######### print only unique alignements


#my %LIN;

#open my $IN2, "<Aligned_TEV_ALL_nolasts";
#open my $OUT2, ">Aligned_TEV_ALL_nolasts_Unique";

#my @ARAW;

#while ( my $line = <$IN2> ) {
#	chomp $line;
#	my @tmp = split (/\t/, $line);
#	my $conc = "$tmp[2]_$tmp[3]";
#	$LIN{"$tmp[2]_$tmp[3]"} = $line;
#}

#foreach my $bla (keys %LIN) {
# print $OUT2 "$LIN{$bla}\n";
#}
