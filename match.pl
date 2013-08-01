#! perl -w 
#find sequence fragments in genome
use strict;
use Bio::SeqIO;
use Bio::AlignIO;
use Bio::Tools::Run::Alignment::Clustalw;

my $genomefile=shift;
my $fragsfile=shift;

my $genome = Bio::SeqIO->new(-file => $genomefile, -format => 'fasta')->next_seq;
my $frags = Bio::SeqIO->new(-file => $fragsfile, -format => 'fasta');
open(OUT, ">seeds.txt");
print OUT "SEQ_NAME\tSTART\tEND\tSTRAND\n";
my $i=1;
while(my $frag=$frags->next_seq){
	print $i."->";
	my ($start,$end)=align($genome,$frag);		
	print $start."..".$end."\n";
	$i++;
}

#make alignment by clustalw
sub align{
	my @aln=@_;
	my @params=('ktuple'=>4,
	    	    'pairgap'=>400,
	    	    'quiet'=>1);
	my $factory=Bio::Tools::Run::Alignment::Clustalw->new(@params);
	my $aln=$factory->align(\@aln);
	my $fragment=$aln[1];
	my $start=$aln->column_from_residue_number($fragment->display_id,1);
	my $end=$aln->column_from_residue_number($fragment->display_id,$fragment->length);
	print $aln->percentage_identity."\t";
	return ($start,$end);
}

