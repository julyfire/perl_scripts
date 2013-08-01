#! perl -w


use strict;
use Bio::SearchIO; 
use Bio::SeqIO;
use Bio::Seq;
use Bio::DB::GenBank;

my $help=<<EOF;

fetch the hited sequences of given blast results

usage: perl $0 parameter_name=parameter_value
parameters:
	i: input file name, blast result format
	e: blast e-value, default 1
	o: output file name, fasta format, default hits.fasta
	f: local blast database file, fasta format, if exists fetch from this file

by weibo, 2011-11-15
EOF

my %parameter;

if(!@ARGV){
	print $help;
	exit;
}

for(@ARGV){
	my($key,$value)=split "=";
	$parameter{$key}=$value;
}

my $infile=$parameter{'i'};
my $e=$parameter{'e'};
$e=1 if !$e;
my $outfile=$parameter{'o'};
$outfile="hits.fasta" if !$outfile;
my $dbfile=$parameter{'f'};

my @hits;

my $in = new Bio::SearchIO(-format => 'blast', 
                           -file   => $infile);
while( my $result = $in->next_result ) {

  while( my $hit = $result->next_hit ) {
    
    while( my $hsp = $hit->next_hsp ) {

	if($hsp->evalue<$e){
		push @hits,$hit->accession;
	}
    }  
  }
}
my $num=scalar(@hits);
print "hit ".$num." sequences\n";

if($dbfile){
	fetchFromLocal($dbfile,@hits);
}
else{
	fetchFromNCBI(@hits);
}

sub fetchFromLocal{
	my $file=shift @_;
	my @acc=@_;
	print "Now start to fetch these sequence from ".$file."\n";
	my %db;
	my $seqio=Bio::SeqIO->new(-file=>$file,-format=>'fasta');
	while(my $seq=$seqio->next_seq){
		$db{$seq->accession_number}=$seq;
	}
	my $outseq=Bio::SeqIO->new(-file=>$outfile,-format=>'fasta');
	my $i=1;
	for(@acc){
		next if !exists($db{$_});		
		$outseq->write_seq($db{$_});
		print $i."/".$num."\n";
		$i++;
	}
	
}

sub fetchFromNCBI{
	my @acc=@_;
	print "Now start to fetch these sequence from NCBI:\n";
	
	my $db=new Bio::DB::GenBank;
	my $outseq=Bio::SeqIO->new(-file=>">hits.fasta",-format=>'fasta');
	my $i=1; 
	for (@acc){
		my $seqio=$db->get_Stream_by_acc($_);
		$outseq->write_seq($seqio->next_seq);
		print $i."/".$num."\n";
		$i++;
	}
}


