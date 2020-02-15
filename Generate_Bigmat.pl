#!/usr/bin/perl -w
use warnings;
use strict;


my @OryMat;

open my $INM, "<$ARGV[0]";
open my $OUT4, ">$ARGV[1]";
while ( my $line = <$INM> ) {chomp $line; if ($line !~ '#') {push @OryMat, $line};}

my @BigScores;
my $max;

for (my $i = 0; $i <= $#OryMat; $i++){
	$BigScores[$i][$i] = 0; 
	my @SS = split (/\s/, $OryMat[$i]);
	my $j = $i + 1;
	$BigScores[$j][$j] = 0;
for (my $k = 0; $k <= $#SS; $k++){
	$BigScores[$i][$j] = $SS[$k];
	$BigScores[$j][$i] = $SS[$k];
	$j++;   ####Some problem here correct!!
	}
	}

print $OUT4 "Reps";

for (my $i = 0; $i <= $#OryMat; $i++){
print $OUT4 "\t$i"
}

print $OUT4 "\n";

for (my $i = 0; $i <= $#OryMat; $i++){
print $OUT4 "$i";
for (my $j = 0; $j <= $#OryMat; $j++){
print $OUT4 "\t$BigScores[$i][$j]"
}
print $OUT4 "\n";
}
