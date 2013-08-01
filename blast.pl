#! /usr/bin/perl -w
#handle blast report file
#by wb 12-02-2009
use strict;
use Bio::SearchIO; 
use Bio::SeqIO;
use Bio::Seq;
use Bio::DB::GenBank;

my $infile=shift;
my $outfile=$infile.".out";




my $in = new Bio::SearchIO(-format => 'blast', 
                          -file   => $infile);
open(OUT,">",$outfile);

my $cc=0;
while( my $result = $in->next_result ) {
 print OUT $result->query_name."\n";
  while( my $hit = $result->next_hit ) {
    my $report="\t".$hit->name.": ";
    while( my $hsp = $hit->next_hsp ) {
      
      if($hsp->evalue<1){
	print OUT $report."(".$hsp->start('hit')."..".$hsp->end('hit').")\tE=".$hsp->evalue."\n";
print "strand: ".$hsp->strand('hit')."\n";
print "start: ".$hsp->start('hit')."\n";
print "end: ".$hsp->end('hit')."\n\n";
      }
    }  
  }
  $cc++;
  print OUT "-----------------\n" if $cc%2==0;
}
exit;
