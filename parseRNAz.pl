#! perl -w
#parse RNAz results

use strict;
use Getopt::Std;

my %opt=();
my $infile="";
my $basename="";
my $startbp=1;
my $start=0;
my $end=0;
my $prob=0.5;
my $direction=0;
my $bpdist=60;
my @steps;

system("clear");
&checkPara();
&load($infile,$startbp);

for my $i (&filter($start,$end,$prob,$direction)){
	print &getStart($i)."\n";
	print &getSeq($i)."\n";
	#print &getStr($i)."\n";
	print &dilute(&getStr($i),$bpdist)."\n";
}

sub checkPara(){
	my $opt_string="b:s:e:p:d:D:h";
	getopts($opt_string,\%opt);
	&usage() if $opt{'h'};
	$infile=shift @ARGV;
	if(!$infile){
		print "Error: need RNAz results file!\n";
		print "Add -h for detail usage\n";
		exit;
	}
	$basename=${[split /\./,$infile]}[0];
	
	
	$startbp=$opt{'b'} if $opt{'b'};
	$start=$opt{'s'} if $opt{'s'};
	$end=$opt{'e'} if $opt{'e'};
	$prob=$opt{'p'} if $opt{'p'};
	$direction=$opt{'d'} if $opt{'d'};
	$bpdist=$opt{'D'} if $opt{'D'};
}
		

sub load(){

	my ($file,$start)=@_;

	open(IN,$file);

	my $i=0;
	$/="############################  RNAz 1.0  ##############################";
	<IN>;
	while(<IN>){
		chomp;
		my $step;
		if($_=~/direction: (.*?)\n.*probability: (.*?)\n.*>.*?(\d+)-(\d+)\n.*>consensus\n(.*?)\n(.*?) \((.*?) /s){
			#print "match $1 $2 $3-$4\n$5\n$6\n$7"."\n";
			$step={"direction"=>$1,
		        	"prob"=>$2,
		        	"start"=>$3+$start,
		        	"end"=>$4+$start,
		        	"seq"=>$5,
		        	"str"=>$6,
		        	"e"=>$7};
			$steps[$i]=$step;
		}
		else{
			print "no match!\n";
		}
		$i++;
	}
	$/="\n";
	close IN;
	print "Load $i fragments!\n";
	
}

sub printProb(){
	my ($bn)=@_;
	open(PROB,">$bn.prob");
	for my $item (@steps){
		print PROB ($item->{'start'})."\t".$item->{'prob'}."\n" if $item->{'direction'} eq "forward";
	}
}

sub filter(){
	my ($start,$end,$prob,$direction)=@_;
	my @hit;
	for my $r (@steps){
		if ($r->{'start'}<=$start and
			$r->{'end'}>=$end and
			$r->{'prob'}>$prob){
			push(@hit,$r) if($r->{'direction'} eq 'forward' and $direction==0);
	 	       push(@hit,$r) if($r->{'direction'} eq 'reverse' and $direction==1);
	 	       push(@hit,$r) if($direction==2);    
	 	}
	}
	return @hit;
	
}

sub getSeg(){
	my ($i)=@_;
	return $steps[$i];
}
sub getSeq(){
	my ($item)=@_;
	return $item->{'seq'};
}

sub getStr(){
	my ($item)=@_;
	return $item->{'str'};
}

sub getE(){
	my ($item)=@_;
	return $item->{'e'};
}

sub getStart(){
	my ($item)=@_;
	return $item->{'start'};
}

sub getEnd(){
	my ($item)=@_;
	return $item->{'end'};
}

sub dilute(){
	
	my ($str,$bpdist)=@_;

	my @loc;
  	my @base=split //,$str;
  
  	my($i,$j,$c,$value,$sum,$lo)=(0,0,0,0,0,0);
  	
  	for($i=0;$i<$#base;$i++){
    		next unless $base[$i] eq '(';
    		for($j=$i+1;$j<=$#base;$j++){
      			next unless $base[$j] eq ')';
      			$sum=0;
      			for($c=$i;$c<=$j;$c++){
        		$value=0;
        		$value=1 if $base[$c] eq '(';
        		$value=0 if $base[$c] eq '.';
        		$value=-1 if $base[$c] eq ')';
        		$sum+=$value;
        		} 
      			next unless $sum==0;
      			if($j-$i>=$bpdist){
        			$base[$i]='.';
        			$base[$j]='.';
        		}
      			else{
        			$loc[$lo++]=[$i,$j];
        			$i=$j;
        		}
      			last;
    		}
  	}
  	$str=join '',@base;
  	return $str;
}
	

sub usage(){
  print <<"EOF";
  
  This program can parse the RNAz results. 
  
  usage: $0 [parameters] results.rnaz
  
  Parameters:
     
    -b number of start of the sequence
    -s number of start of the segment of interest one specifies
    -e number of end of the segment of interest one specifies
    -p the P value that one specifies (0.5)
    -d the direction of the sequence (0), 0:forward, 1:reverse, 2:two direction
    -D base pairing distance (60)
    -h help
  
  example: $0 [parameters] seq.fasta
  
EOF

exit;
}	
