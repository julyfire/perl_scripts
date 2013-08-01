#! perl -w
use strict;
use Bio::Seq;
use Bio::Tools::Run::RemoteBlast;
my $prog = 'blastn';
my $db   = 'nr';
my $e_val= '1e-10';

my @params = ( '-prog' => $prog,
         '-data' => $db,
         '-expect' => $e_val,
         '-readmethod' => 'SearchIO' );

my $factory = Bio::Tools::Run::RemoteBlast->new(@params);

  #change a paramter
$Bio::Tools::Run::RemoteBlast::HEADER{'ENTREZ_QUERY'} = 'Dengue virus 1';

  #remove a parameter
#delete $Bio::Tools::Run::RemoteBlast::HEADER{'FILTER'};

my $v = 1;
  #$v is just to turn on and off the messages

my $seq=Bio::Seq->new(-display_id=>"D00502",
			-seq=>"TCGAGATGTCCAACACAAGGAGAAGCCACGCTGGTGGAAGAACAGGACACGAACTTTGTGTGCCGACGAACGTTCGTGGACAGAGGCTGGGGCAATGGTTGTGGGCTATTCGGAAAAGGTAGCTTAATAACGTGTGCTAAGTTTAAGTGTGTGACAAAACTGGAAGGAAAGATAGTC",
			-desc=>"test",
			-alphabet=>"dna");


    #Blast a sequence against a database:
my $r = $factory->submit_blast($seq);
    #my $r = $factory->submit_blast('amino.fa');

print STDERR "waiting..." if( $v > 0 );
while ( my @rids = $factory->each_rid ) {
	foreach my $rid ( @rids ) {
        	my $rc = $factory->retrieve_blast($rid); #Bio::SearchIO object
        	if( !ref($rc) ) { #if still searching
          		if( $rc < 0 ) { #if there's no hits
            			$factory->remove_rid($rid);
          		}
          		print STDERR "." if ( $v > 0 );
          		sleep 5;
        	} else {
          		my $result = $rc->next_result();

          		#save the output
          		my $filename = $result->query_name()."\.out";
          		$factory->save_output($filename);
          		$factory->remove_rid($rid);

          		print "\nQuery Name: ", $result->query_name(), "\n";
          		while ( my $hit = $result->next_hit ) {
            			next unless ( $v > 0);
            			print "\thit name is ", $hit->name, "\n";
				print "\thit length is ", $hit->length, "\n";
            			while( my $hsp = $hit->next_hsp ) {
              				print "\t\tscore is ", $hsp->score, "\n";
					#print "strand: ".$hsp->strand('hit')."\n";
					#print "start: ".$hsp->start('hit')."\n";
					#print "end: ".$hsp->end('hit')."\n\n";
					#print "length: ".$hsp->length,"\n";
					#print "hit: ".$hsp->seq('hit')->seq,"\n";
					#print "hit: ".$hsp->seq_str('hit'),"\n";
					#print "query: ".$hsp->query_string,"\n";
					print "hit: ".$hsp->hit_string,"\n";
            			}
          		}
        	}
	}
}

