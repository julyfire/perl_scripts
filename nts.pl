#! perl -w
# the program can calculate the number of given oligonucleotides
# by weibo 2010-12-7
use strict;

my $infile="";
my @ntlist;
my $outfile="";
&checkcmd();

my %seqs=&readfasta($infile);

open(OUT, ">",$outfile);

#print the column names
print OUT "sequence\t";
for my $nt (@ntlist){
	print OUT $nt."\t";
}
print OUT "\n";

for my $seq (keys %seqs){
	print OUT $seq."\t";
	for my $nt (@ntlist){
		my $num=&content($seqs{$seq},$nt);
		print OUT $num."\t";
	}
	print OUT "\n";
}

sub content{
	my($seq,$nt)=@_;
	my $length=length($seq);
	my @nts=&cut($seq,length($nt));
	my $count=0;
	for my $base (@nts){
		$count++ if $base eq $nt;
	}
	return $count;
}

sub cut{
	my($seq,$window)=@_;
	my @unit;
	for(my $i=0;$i<=length($seq)-$window;$i++){
		$unit[$i]=substr($seq,$i,$window);
	}
	return @unit;
}
	
sub readfasta{
	my($fasta)=@_;
	my %seqs;
	open(IN, $fasta) or die "Error: cannot open the file ".$fasta."\n";
	$/=">";
	<IN>;
	while(<IN>){
		chomp;
		s/\r//gs;
		my $s=index($_,"\n");
		my $id=substr($_,0,$s);
		my $seq=substr($_,$s);
		$seq=~s/[\r\n\s\-]+//gs;
		$seqs{$id}=uc $seq;
	}
	$/="\n";
	close IN;
	return %seqs;
}	

sub checkcmd{
	if(scalar @ARGV==2){
		$infile=$ARGV[0];
		@ntlist=split /,/, uc $ARGV[1];
		$outfile=${[split /\./,$infile]}[0].".xls";
	}
	elsif(scalar @ARGV==3){
		$infile=$ARGV[0];
		@ntlist=split /,/, uc $ARGV[1];
		$outfile=$ARGV[2];
	}
	else{
		print "The program can calculate the number of given oligonucleotides in fasta format sequences\n"; 
		print "\nUSAGE: perl $0 fasta_file oligonucleotides [output_filename]\n\n";
		print "e.g. input:\n";
		print "            perl $0 seq.fasta gc,at\n\n";
		print "     will calculate the number of gc and at in file seq.fasta,\n";
		print "     the result will be saved in file seq.xls\n\n";
		print "one can specify the third parameter as the output file name\n";
		print "more than one oligonucleotides should be separated by comma\n\n";
		exit 0;
	}
	
}	
	
	
