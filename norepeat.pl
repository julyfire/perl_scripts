#! perl -w
#remove repeat sequence from fasta MSA
#usage: perl norepeat.pl infile

use strict;

my $infile=shift;
open(IN,$infile);

my %seqs;

$/=">";
while(<IN>){
	chomp;
	next if $_ eq "";
	my($id,$seq)=/(.*?)\n(.*)/s;	
	$seqs{$seq}=$id;
}
$/="\n";
close IN;

$infile=~s/(.*)\..*/$1/;
my $outfile=$infile."-norepeat.fasta";
open(OUT,">".$outfile);
for my $seq (keys %seqs){
	print OUT ">".$seqs{$seq}."\n".$seq;
}
close OUT;
