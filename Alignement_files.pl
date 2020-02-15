#!/usr/bin/perl -w
use warnings;
use strict;

#usage perl Alignement_files.pl TEV2_ALL_Repeatscode_withstars.txt Multiple_align

if ( -d "ALIGNS") {} else { system ("mkdir ALIGNS")};
if ( -d "ALIGNS/ForFUNC") {} else { system ("mkdir ALIGNS/ForFUNC")};

#Store matriX in hashes

my $files = "BigRepDist.mat";
my %DIST; #Hash of hashes, each hash being a multidimendional hassh with two values

open my $MAT, "<$files" or die "Can not open $files, file doesn't exist\n";
print "reading $files\n";
my $A = 0;
while ( my $line = <$MAT> ) {
 my $B = 0;
	chomp $line;
	if ($line =~ "Reps") {} else {
	my @ARR = split (/\t/, $line);
	for (my $j = 1; $j <= $#ARR; $j++) {
	$DIST{$A}{$B} = $ARR[$j]; 
	#print "$A vs $B\n";
	$B++;
 	} #each line
	$A++;

	} #all lines except first
} #inside each file
close $MAT;



#Store RVDs
my %RVDs;

open my $IN0, "<$ARGV[0]";
my @FRAW;

while ( my $line = <$IN0> ) {
	chomp $line;
	push @FRAW, $line;
}

for (my $i =0; $i <= $#FRAW; $i++) {#mod ALPQ
		my @ARR = split (/\t/, $FRAW[$i]);
		my $R = substr ($ARR[1], 11, 2);
		$RVDs{$ARR[0]} = $R;
}


##############READ multiple alignments and store them in hash

open my $IN1, "<$ARGV[1]";
my @ARAW;

while ( my $line = <$IN1> ) {
	chomp $line;
	push @ARAW, $line;
}

my %GAL;
my $G;
for (my $i =0; $i <= $#ARAW; $i++) {#mod ALPQ
		if ($ARAW[$i] =~ "Group") {$G = (split (/\s/, $ARAW[$i]))[2];}
		else {
		my @tmp = split /\s+/, $ARAW[$i];
		my $s = join ("\t", @tmp);
		push @{ $GAL{$G} }, $s;
		}
}



foreach my $Group (sort keys %GAL) {
			my @Sarr =  @{ $GAL{$Group} };
			open my $OUT0, ">ALIGNS/$Group.Reps";
			open my $OUT1, ">ALIGNS/$Group.RVDs";
			open my $OUT2, ">ALIGNS/$Group.Dists";
			open my $OUT3, ">ALIGNS/ForFUNC/$Group.RVDs";
			my $first = "$Sarr[0]";
				my  @FRep = split (/\t/, $first); my @FRep2; my @FRep3; #reference arrays for distances
				if (defined $Sarr[1]) {@FRep2 = split (/\t/, $Sarr[1])};
				if (defined $Sarr[2]) {@FRep3 = split (/\t/, $Sarr[2])};
			for (my $i =0; $i <= $#Sarr; $i++) {
			print $OUT0 "$Sarr[$i]\n";
			my  @Reps = split (/\t/, $Sarr[$i]);
			print $OUT1 "$Reps[0]";print $OUT2 "$Reps[0]";print $OUT3 "$Reps[0]\t";
			for (my $j =1; $j <= $#Reps; $j++) {
			if ($Reps[$j] eq "-") {print $OUT1 "\t-"; print $OUT2 "\tNA"; }
			else {print $OUT1 "\t$RVDs{$Reps[$j]}"; print $OUT3 "$RVDs{$Reps[$j]}-";
					if (defined $DIST{$Reps[$j]}{$FRep[$j]}) {
					  print $OUT2 "\t$DIST{$Reps[$j]}{$FRep[$j]}"
															} else {print $OUT2 "\t$DIST{$Reps[$j]}{$FRep2[$j]}"}
				
				}
			
											}
			print $OUT1 "\n";print $OUT2 "\n";print $OUT3 "\n";
}
}


