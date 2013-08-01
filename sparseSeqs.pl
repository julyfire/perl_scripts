#! perl -w
# remove sequences with pairwise identity greater than the given cutoff
# usage: perl sparseSeqs.pl inputfile cutoff
# by weibo 2013-7-22

use strict;
use Bio::SeqIO;
use Bio::AlignIO;

my $infile=shift;
my $cutoff=shift;
my $outfile=${[split /\./,$infile]}[0]."_$cutoff.fasta";

my $in=Bio::AlignIO->new(-file=>$infile,-format=>"fasta");
my $seqout=Bio::SeqIO->new(-file=>">".$outfile,-format=>"fasta");

print "Input file: ".$infile."\n";
print "identity cutoff: ".$cutoff."%\n";

while(my $aln=$in->next_aln){
	my $n=$aln->num_sequences;
	my @left=(1..$n);
	my @tmp=();
	for(my $i=0;$i<$#left;$i++){
		for(my $k=0;$k<=$i;$k++){
			push @tmp, $left[$k];
		}
		for(my $j=$i+1;$j<=$#left;$j++){
			my $pid=$aln->select_noncont($left[$i],$left[$j])->percentage_identity;	
			if($pid<$cutoff){
				push @tmp,$left[$j];								
			}
			else{
				my $s1=$aln->get_seq_by_pos($left[$i])->id;
				my $s2=$aln->get_seq_by_pos($left[$j])->id;
				print "remove ".$s2." having ".$pid."% of identity with ".$s1."\n";
			}
		}
		@left=@tmp;
		@tmp=();
	}
	print "keep ".scalar(@left)." sequences\n";
	for my $i (@left){
		$seqout->write_seq($aln->get_seq_by_pos($i));
	}
	print "see output file in ".$outfile."\n";
}
