#!perl -w
#extracts specified region from genbank file
#output: .fasta file.
#by wb 2010-09-14

use strict;
use Bio::SeqIO;
use Bio::Seq;

my $infile=shift;
my $region=shift;
if (!$infile){
	usage();
	exit;
}

my $seqio=Bio::SeqIO->new(-file=>$infile, -format=>'genbank');
my $outfile=${[split /\./, $infile]}[0].".fasta";
my $seqout=Bio::SeqIO->new(-file=>'>'.$outfile, -format=>'fasta');
my $new_seq=Bio::Seq->new();

my ($acc,$organism,$country,$date);
while(my $seq=$seqio->next_seq){
	my($start,$end)=(0,0);
	($acc,$organism,$country,$date)=('','','','');
	$acc=$seq->accession_number;
	print $acc."\n";
	for my $feat ($seq->get_SeqFeatures) {
		if ($feat->primary_tag eq 'source'){
			$organism=join '',$feat->get_tag_values('organism') if ($feat->has_tag('organism'));
			$country=join '',$feat->get_tag_values('country') if ($feat->has_tag('country'));
			$date=join '',$feat->get_tag_values('collection_date') if ($feat->has_tag('collection_date'));
		}
		if($feat->primary_tag eq 'Protein' and $feat->has_tag('product') and (join '',$feat->get_tag_values('product'))=~/NS5/i){
			$start=$feat->location->start;
			$end=$feat->location->end;			
		}
		elsif($feat->primary_tag eq 'mat_peptide' and $feat->has_tag('product') and (join '',$feat->get_tag_values('product'))=~/NS5/i){
			$start=$feat->location->start;
			$end=$feat->location->end;
		}
		elsif($feat->primary_tag eq 'Region'){
			if($feat->has_tag('region_name') and (join '',$feat->get_tag_values('region_name'))=~/MTases/i){
				$start=$feat->location->start;
				$end=$feat->location->end;
			}
			if($feat->has_tag('region_name') and (join '',$feat->get_tag_values('region_name'))=~/NS5/i){
				$start=$feat->location->start if $start==0;
				$end=$feat->location->end;
			}	
		}
	}
	$new_seq->display_id($acc);
	$new_seq->desc($organism."_".$country."_".$date);
	$new_seq->seq($seq->subseq($start,$end));

	$seqout->write_seq($new_seq);
}

sub usage{
	print "\nextracts specified region from genbank file and save as *.fasta file.\n";
	print "usage: perl $0 input_file\n";
	print "--input_file: genbank format file (.gb)\n";
}
