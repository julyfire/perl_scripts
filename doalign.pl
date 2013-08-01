 #! perl -w
 use strict;
 use Bio::SeqIO;
 
 my @a=<./primers/*>;
 for(my $i=0;$i<scalar(@a);$i++){
    $i+=1;
    my $seqp=Bio::SeqIO->new(-file=>"./primers/".$i.".fasta",-format=>'fasta');
    my $primer=$seqp->next_seq;
    my $pid=$primer->id;
    
    my $dir="./aligns/";
    
    open(OUT,">".$dir.$i.".fasta");
    print OUT ">".$primer->id."\n";
    print OUT $primer->seq."\n";
    
    system("needle -asequence ./primers/".$i.".fasta -snucleotide1 1 -sformat1 fasta -bsequence consensus.fasta -snucleotide2 1 -sformat2 fasta -outfile ".$dir.$i.".needle -aformat3 pair -gapopen 10.0 -gapextend 0.5");
    
    my $input=$dir.$i.".needle";
    use Bio::AlignIO;
    # read in an alignment from the EMBOSS program water
    my $in = Bio::AlignIO->new(-format => 'emboss',
                              -file   => $input);
    while( my $aln = $in->next_aln ) {
        my $primer=$aln->get_seq_by_pos(1); 
        my $start=$aln->column_from_residue_number( $primer->id, 2 )-1;
        my $end=$aln->column_from_residue_number($primer->id,$primer->end-1)+1;
        my $aln2=$aln->slice($start,$end);
        
        
        my $temp=$aln2->get_seq_by_pos(2);
        print OUT ">".$temp->id."\n";
        print OUT $temp->seq."\n";
    
    }
    close OUT;
    unlink $input;
    system("edialign  -sequences ".$dir.$i.".fasta  -sformat1 fasta -snucleotide1 -nucmode n -norevcomp  -overlapw 1 -linkage UPGMA -maxfragl 40 -nofragmat  -fragsim 4 -noitscore  -threshold 0.0 -nomask  -nodostars  -starnum 4 -osformat fasta -odirectory2 ".$dir." -osname3 ".$i." -osdirectory3 ".$dir." -auto");
    unlink $dir.$i.".edialign";

    $i-=1;
 }
 
