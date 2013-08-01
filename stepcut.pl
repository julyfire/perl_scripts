#! perl -w
#cut a part of an fasta alignment by specified window and step
#weibo 2010-7-12
use strict;
use Bio::SeqIO;
use Bio::Seq;
#use Bio::AlignIO;
use Getopt::Std;

system("clear");

my ($infile,$bn,$start,$end,$step,$win,$length,$dir);
my %opt=();

###-------------start check parameters-------
my $opt_str="i:b:e:s:w:d:h";
getopts($opt_str,\%opt);
&usage() if $opt{h};
if($opt{i}){
	$infile=$opt{i};
	$bn=${[split /\./,$infile]}[0];
	$length=Bio::SeqIO->new(-file=>$infile)->next_seq->length;
}else{
	print "E: must assign a fasta alignment file!\n";
	print "see help by -h\n";
	exit;
}
if($opt{b}){
	$start=$opt{b};
}else{
	print "W: start point is not assigned, using default 1\n";
	$start=1;
}
if($opt{e}){
	$end=$opt{e};
}else{
	print "W: stop point is not assigned, using default sequence end\n";
	$end=$length;
}
if($opt{s}){
	$step=$opt{b};
}else{
	print "W: step size is not assigned, using default 40\n";
	$step=40;
}
if($opt{w}){
	$win=$opt{w};
}else{
	print "W: window size is not assigned, using default 200\n";
	$win=200;
}
if($opt{d}){
	$dir=$opt{d};
	system("mkdir $dir");
}else{
	$dir=".";
}
###-------------end of parameters checking----------------

my $n=int(($end-$start+1-$win)/$step+1);

my $seqin=Bio::SeqIO->new(-file=>$infile);
my @pos;
my @seqout;
for(my $i=0;$i<$n;$i++){
	$pos[$i]=[$start+$step*$i,$start+$step*$i+$win-1];
	$seqout[$i]=Bio::SeqIO->new(-file=>">$dir/".$bn."_".$pos[$i][0]."-".$pos[$i][1].".fasta",
					-format=>"fasta");
}

while(my $seq=$seqin->next_seq){
	for(my $i=0;$i<$n;$i++){
		$seqout[$i]->write_seq($seq->trunc($pos[$i][0],$pos[$i][1]));
	}
}

sub usage(){
  print <<"EOF";
  
  This program cuts a part of an fasta alignment by specified window and step.
  by weibo 2010-09-15 
  
  Parameters:
    -i the input file (.fasta)
    -b the start point of cutting, default 1
    -e the stop point of cutting, default the length of the sequence
    -s step size, default 40
    -w window size, default 200
    -d the directory to save the out file
    -h help
  
  example: $0 -i seq.fasta
  
EOF

exit;
}
