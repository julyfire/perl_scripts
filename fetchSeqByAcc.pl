#!usr/bin/perl -w
# by WB 11.25.08
# usage: >perl acc.pl infile outfile format
use strict;
use Bio::SeqIO;
use Bio::DB::GenBank;
my $use="This program is used to retrive ncbi file with accession number.\nUsage:\n      perl acc.pl infile_name outfile_name outfile_format\n";
my $infile=shift or die print $use;
my $outfile=shift or die print $use;
my $format=shift or die print $use;
open(INFILE,$infile)or die "cannot open the file: $infile\n";
my @acc=<INFILE>;
my $i=0;
for (@acc){
my $db=new Bio::DB::GenBank;
my $seqio=$db->get_Stream_by_acc($_);
my $outseq=Bio::SeqIO->new(-file=>">>$outfile",-format=>$format);
while(my $seq=$seqio->next_seq){
$outseq->write_seq($seq);
}
$i++;
print "fetched $i\n";
}
print "\n$i entries are fetched!\n";
exit;
