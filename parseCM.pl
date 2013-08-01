#! perl -w
#cut an alignment fasta file into fragments(default 400bp) with breakpoints as the mid-base, then feed CMfinder to find some secondary structure pattern, at last draw the structure through VARNA.
#by weibo 2010-09-17 

use strict;
use Bio::SeqIO;
use Bio::Seq;

my $infile=shift;
my $bps=shift;

if(!$infile or !$bps){
	print "usage: $0 infile breakpoint1,breakpoint2,...\n";
	exit;
}
my @breakpoint=split /,/,$bps;
my $extend=200;
my $bn=substr($infile,0,rindex($infile,'.'));
my $dir="./temp";
my $seqStart=4608;
system("mkdir $dir");

for my $bp (@breakpoint){
	my($start,$end)=($bp-$extend,$bp+$extend);
	my $outname=$dir."/".$bn."_".$start."-".$end;
	
	# cut alignment
	my $seqin=Bio::SeqIO->new(-file=>$infile);
	my $seqout=Bio::SeqIO->new(-file=>">$outname.fasta", -format=>"fasta");
	my $segStart=$start-$seqStart+1;
	my $segEnd=$end-$seqStart+1;
	while(my $seq=$seqin->next_seq){
       	$seqout->write_seq($seq->trunc($segStart,$segEnd));
        }
        
       # cmfinder
       my $result=`cmfinder.pl $infile`;
	open(OUT,">$outname.summary");
	print OUT $result;
	close OUT;
	
	my @rna=($result=~/Seq_.*?_(\d+)_(\d+).+?\n(.*?)\n(.*?)\n/gs);
	
	# VARNA
	for(my $i=0;$i<=$#rna;$i+=4){
		my $strStart=$rna[$i]+$start-1;
		my $strEnd=$rna[$i+1]+$start-1;
		my $seq=$rna[$i+2];
		my $str=$rna[$i+3];
		my $strOut=$bn."_".$strStart."-".$strEnd.".eps";
		system("java -cp /home/weibo/Desktop/VARNA fr.orsay.lri.varna.applications.VARNAcmd -sequenceDBN '$seq' -structureDBN '$str' -o $strOut -flat true -title '$bn' -startNum $strStart");
	} 
}



