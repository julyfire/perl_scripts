#!perl -w
#extracts some features from genbank file
#output: .txt file including accession,organism,country and date in each line.
#by wb 2010-09-14

use strict;
use Bio::SeqIO;
use Bio::Seq;
use Bio::DB::GenBank;

system('clear');

my $accfile=shift;
if (!$accfile){
		print "usage: perl $0 accfile\n";
		exit;
	}
open(ACCFILE, $accfile);
my @acc=<ACCFILE>;
my $i=1;
for my $acc (@acc){
	chomp($acc);
	$acc=~s/\r//;
	next if $acc eq "";
	getGb($acc);
	getCDS($acc);
	print "$i:$acc\n";
	$i++;
}

####################sub routine#########################
sub getGb{
	my($acc)=@_;
	my $db=new Bio::DB::GenBank;
	my $seqio=$db->get_Stream_by_acc($acc);
	my $outseq=Bio::SeqIO->new(-file=>">$acc.gb",-format=>"genbank");
	while(my $seq=$seqio->next_seq){
		$outseq->write_seq($seq);
	}
}

sub getCDS{
	my ($acc)=@_;
	my $infile=$acc.".gb";
	#if (!$infile){
	#	print "usage: perl $0 seq.gb\n";
	#	exit;
	#}

	my $seqio=Bio::SeqIO->new(-file=>$infile, -format=>'genbank');
	my $outfile="out_CDS.fasta";
	my $outfile2="out_protein.fasta";
	open(OUT, ">>$outfile");
	open(OUT2, ">>$outfile");
	my ($organism,$country,$date);
	while(my $seq=$seqio->next_seq){
		($organism,$country,$date)=('','','');
		#$acc=$seq->accession_number;
		
		for my $feat ($seq->get_SeqFeatures) {
			if ($feat->primary_tag eq 'source'){
				$organism=${[$feat->get_tag_values('organism')]}[0] if ($feat->has_tag('organism'));
				$country=${[$feat->get_tag_values('country')]}[0] if ($feat->has_tag('country'));
				$date=${[$feat->get_tag_values('collection_date')]}[0] if ($feat->has_tag('collection_date'));
		
				print OUT ">$acc"."_".$organism."_".$country."_".$date."\n";
				print OUT2 ">$acc"."_".$organism."_".$country."_".$date."\n";
			}
			if($feat->primary_tag eq 'CDS'){
				print OUT $feat->seq->seq."\n";
				print OUT2 ${[$feat->get_tag_values('translation')]}[0]."\n" if ($feat->has_tag('translation'));
			}
		}
	}
	system("rm $acc.gb");
}


