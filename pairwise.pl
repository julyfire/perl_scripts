#! perl -w 
use strict;
use Bio::SeqIO;
use Bio::AlignIO;




my $seq1 = Bio::SeqIO->new(-file => $ARGV[0], -format => 'fasta')->next_seq;

my $seq2 = Bio::SeqIO->new(-file => $ARGV[1], -format => 'fasta')->next_seq;

my @aln=($seq1,$seq2);
my $aln;
my $factory;

=head1
#make alignment by bioperl-ext package
use Bio::Tools::dpAlign;
$factory = new dpAlign(-match => 3,
                     -mismatch => -1,
                     -gap => 3,
                     -ext => 1,
                     -alg => Bio::Tools::dpAlign::DPALIGN_ENDSFREE_MILLER_MYERS);

my $out = $factory->pairwise_alignment($seq1->next_seq, $seq2->next_seq);
my $alnout = Bio::AlignIO->new(-format => 'pfam', -fh => \*STDOUT);
$alnout->write_aln($out);
=cut


#make alignment by blast2
use Bio::Tools::Run::StandAloneBlast;
$factory = Bio::Tools::Run::StandAloneBlast->new(-program => 'blastn',
                                                   -outfile => 'bl2seq.out');
my $bl2seq_report = $factory->bl2seq($seq1, $seq2);
# Use AlignIO.pm to create a SimpleAlign object from the bl2seq report
$aln = Bio::AlignIO->new(-file=> 'bl2seq.out',-format => 'bl2seq')->next_aln();


=head
#make alignment by clustalw
use Bio::Tools::Run::Alignment::Clustalw;
my @params=('ktuple'=>4,
	    'pairgap'=>400,
	    'quiet'=>1);
$factory=Bio::Tools::Run::Alignment::Clustalw->new(@params);
$aln=$factory->align(\@aln);
my $alnout = Bio::AlignIO->new(-format => 'fasta', -fh => \*STDOUT);
$alnout->write_aln($aln);
=cut

