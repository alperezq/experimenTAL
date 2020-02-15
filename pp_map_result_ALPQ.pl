
#!/usr/bin/perl
$argc = @ARGV;


if ($argc < 2)   
{
    print "usage: \n perl  pp_map_result.pl MSATcompareOutfile SequenceFile\n";
    exit;
}
my $infile = $ARGV[0];
my $seq_file = $ARGV[1];
#my $output_tex = $ARGV[2];

my %ID;
my $count = 0;

open(IN, "<$infile");
my @lines=<IN>;
close(IN);

open(seqIN, "<$seq_file");
my @preSeqs=<seqIN>;
close(seqIN);
my @Seqs;
foreach $Seq (@preSeqs){
    if($Seq=~/>/){
	chomp $Seq;
	$Seq =~ s/>//g;
	$ID{$count} = $Seq;
	$count++; 
    }
    else {push (@Seqs,$Seq)};
}
foreach $Seq (@Seqs){
    if($Seq=~/>/){	
	next;
    }
}

#exit; 
my $stringS="";
my $stringR="";





### part 1 read alignment
my $flag=0;

#my $counterS=0;
#my $counterR=0;
my $line_beg=0;
my $line_end=0;
my $comp_counter=0;

$line_n=0;
foreach $line2 (@lines){
    if($line2=~/operations/){
	$flag=1;
	$line_beg=$line_n;
	#next;
    }
    if($line2=~/Score of aligning/){
	$line2=~/Seq:(\d+).+Seq:(\d+)/;
	$seq1_no=$1;
	$seq2_no=$2;
	
	$flag=1;
	$line_end=$line_n;	
	$stringS=$Seqs[$seq1_no];
	$stringR=$Seqs[$seq2_no];

	$comp_counter=$seq1_no."x".$seq2_no;
	#print "$comp_counter\n";
	#print "$ID{$seq1_no} vs $ID{$seq2_no}\n" ;
	my @AL = draw_alignment($line_beg,$line_end,$comp_counter);
	print "$ID{$seq1_no}\t$ID{$seq2_no}\t$AL[0]\t$AL[1]\n";
	$comp_counter++;
	#next;
    }
    
    $line_n++;
}

exit;
####

sub draw_alignment{
    $begin=$_[0];
    $end=$_[1];
    $outputfilecounter=$_[2];

    my @seqS=split(" ",$stringS);
    my @seqR=split(" ",$stringR);	

    my $AliS;
    my $AliR;
    my $operation;
    my $line;
    for ($i=$begin;$i<=$end;$i++){
	
	$line=$lines[$i];
	if($line=~/operations/){
	    $flag=1;
	    next;
	}
	if(length($line)<3){
	    next;
	}
	if($flag==1){
	    if($line=~/Match.+\((\d+),\s(\d+)\)/){
		if($1==0) {next;}
		$AliS=$AliS." ".$seqS[$1-1]; #modALPQ removealignement start
		$AliR=$AliR." ".$seqR[$2-1];
	    }
	    if($line=~/dup.+S.+\[(\d+)..(\d+)\]/){
		$lb=$1;$rb=$2;
		if($line=~/Left/){$lb=$lb+1;}
		if($line=~/Right/){$rb=$rb-1;}
		for($k=$rb;$k>=$lb;$k--){
		    $AliS=$AliS." ".$seqS[$k-1];
		    $AliR=$AliR." -";
		}
		#$AliS=$AliS." ".$seqS[$1-1];
		#$AliR=$AliR." ".$seqR[$2-1];
	    }
	    if($line=~/dup.+R.+\[(\d+)..(\d+)\]/){
		$lb=$1;$rb=$2;
		if($line=~/Left/){$lb=$lb+1;}
		if($line=~/Right/){$rb=$rb-1;}
		for($k=$rb;$k>=$lb;$k--){
		    $AliR=$AliR." ".$seqR[$k-1];
		    $AliS=$AliS." -";
		}
		#$AliS=$AliS." ".$seqS[$1-1];
		#$AliR=$AliR." ".$seqR[$2-1];
	    }
	}
	$line_no++;
    }
    
    
    @tmp=split(/ /,$AliS);
    @tmp=reverse @tmp;
    
    $AliS2="";
    for($k=0;$k<scalar @tmp;$k++){
	$AliS2=$AliS2." ".$tmp[$k];	
    }

    @tmp=split(/ /,$AliR);
    @tmp=reverse @tmp;
    
    $AliR2="";
    for($k=0;$k<scalar @tmp;$k++){
	$AliR2=$AliR2." ".$tmp[$k];	
    }
#    $AliS2=reverse $AliS;
#    $AliR2=reverse $AliR;

#    $AliS2="\\\$ ".$AliS2;
#    $AliR2="\\\$ ".$AliR2;
    
    #print $AliS2."\n";
    #print $AliR2."\n";
	my @Arr = ("$AliS2", "$AliR2");
	return @Arr;
    

    
}
