#! /usr/bin/perl -w

############ weibo 2009-9-9 14:39 #########################
######### mi & bootstrap ###########################3
use strict;
use Bio::SeqIO;
my $input=shift;
my $outMI=$input."."."MI3";#mi matrix outfile
my $outMIp=$input."."."MIp";#
my $outbs=$input."."."bs";#mi cumulative probibality distribution outfile

my $cutoff=0.0002;

my $seqio=Bio::SeqIO->new(-file=>$input,-format=>'fasta');
open(OUT,">$outMI");
open(OUTO,">$outMIp");
open(OUTBS,">$outbs");
my $length=$seqio->next_seq->length;#alignment length

my $key;
my $value;

my $p;

my $count=0;#seqs number

####################################     figure out the entropy
my @entropy;#store the entropy of each site
my @seqs;#store the seqs matrix,each element is a base

print "start to calculate Entropy!\n";
for(my $i=0;$i<$length;$i++){
    my %sites;
    
    my $sequence;
    $count=0;
    $seqio=Bio::SeqIO->new(-file=>$input,-format=>'fasta');
    while($sequence=$seqio->next_seq){
        my $seq=$sequence->seq;
        $p=substr($seq,$i,1);
	
	$seqs[$count][$i]=$p;
	
        if(exists $sites{$p}){
            $sites{$p}+=1;
        }
        else{
            $sites{$p}=1;
        }
        $count++;
    }
    
    while(($key,$value)=each %sites){
	$entropy[$i]+=-$value/$count*log($value/$count)/log(20);  #...H(a)
	#print "Entropy: ".$key."=>".$value."\n";
        
    }
   
}
print "The entropy calculation has finished!\n";

print "Now start to calculate MI(a,b)\n";
my @mi=&MI(@seqs);
print "MI calculation has finished!\n";

#print out MI matrix
for(my $i=0;$i<$length;$i++){
    for(my $j=0;$j<$length;$j++){
        printf OUT "%2.10f\t",$mi[$i][$j];   
    }
    print OUT "\n";
   
}

=head
##################3 bootstrap ###########################3
print "Input bootstrap number: ";
chomp(my $number=<STDIN>);
my %vmi=&bootstrap($number);
foreach my $key(sort{$a<=>$b} keys %vmi){
    print OUTBS $key."\t".$vmi{$key}."\n";
}

#filter MI matrix by bootstrap probability
for(my $i=0;$i<$length;$i++){
    for(my $j=0;$j<$length;$j++){
	$mi[$i][$j]=filter($mi[$i][$j],$cutoff);
        printf OUTO "%2.10f\t",$mi[$i][$j];   
    }
    print OUTO "\n";
   
}

sub filter{
	my ($v,$cutoff)=@_;
	my @keys=sort{$a<=>$b} keys %vmi;
	
	for my $key(@keys){
		if($key>=$v){
			if(1-$vmi{$key}<$cutoff){
				$v=$vmi{$key};
			}
			else{
				$v=0;
			}
			last;
		}
	}
	return $v;
}
=cut

##################3 bootstrap ###########################3
print "Input bootstrap number: ";
chomp(my $number=<STDIN>);
my @pm=&bootstrap($number);

#filter MI matrix by bootstrap probability
for(my $i=0;$i<$length;$i++){
    for(my $j=0;$j<$length;$j++){
	printf OUTBS "%2.10f\t",$pm[$i][$j];
	if ($pm[$i][$j]<$cutoff){
	        printf OUTO "%2.10f\t",1-$pm[$i][$j];
	}
	else{
		print OUTO "0\t";   
	}
    }
    print OUTO "\n";
    print OUTBS "\n";
}



sub MI{
    my @seqs=@_;
    my @mi;
    #my $count=scalar(@seqs);
    #my $length=$#{$seqs[0]}+1;

   #MI matrix initialization
    for(my $i=0;$i<$length;$i++){
	for(my $j=0;$j<$length;$j++){
	    $mi[$i][$j]=0;
	}   
    }

    for(my $i=0;$i<$length-1;$i++){
	for(my $j=$i+1;$j<$length;$j++){
	    my %pairs;
	
	    for(my $k=0;$k<$count;$k++){
		my $p=$seqs[$k][$i]." ".$seqs[$k][$j];
		if(exists $pairs{$p}){
		    $pairs{$p}+=1;
		}
		else{
		    $pairs{$p}=1;
		}
            
	    }
	    my $hab=0;
	    while(($key,$value)=each %pairs){                      #...H(a,b)
		$hab+=-$value/$count*log($value/$count)/log(20);
	    }
        
	    $mi[$i][$j]=$entropy[$i]+$entropy[$j]-$hab;            #...MI(a,b)
	}
    }
    return @mi;
}
=head
#parameter: number of bootstrap
sub bootstrap{
    my %vmi;#store all mi values of all the bootstraps
    my $num=0;#the total number of mi values of all the bootstraps

    
    srand (time ^ $$ ^ unpack "%L*", `ps axww | gzip -f`);
    
    for(my $n=0;$n<$_[0];$n++){#each bootstrap
    
        my @bseq;
        for(my $j=0;$j<$length;$j++){#each site
            my %rands;
            for(my $i=0;;){
                last if($i==$count);
                my $rn=int(rand($count));#get a seq randomly
                if(exists $rands{$rn}){#if this seq has been used, choose another seq
                    $rands{$rn}=0;
                    next;
                }else{
                    $rands{$rn}=1;
                }
                $bseq[$i++][$j]=$seqs[$rn][$j];#shuffle the bases in one column
                #print $rn." ";
            }
            #print "\n";
        }
    
        my @bmi=&MI(@bseq);#calculate MI of this bootstrap sample
    
	#calculate the times of each mi value
        for(my $i=0;$i<$length-1;$i++){
            for(my $j=$i+1;$j<$length;$j++){
                my $mii=$bmi[$i][$j];
                if(exists $vmi{$mii}){
                    $vmi{$mii}+=1;
                }else{
                    $vmi{$mii}=1;
                }
                $num++;
            }
        }
        print "$n of the ".$_[0]."\n";
    }
    
    #calculate the distribution of probibality of mi value
    my $cum=0;
    foreach my $key(sort{$a<=>$b} keys %vmi){
        $vmi{$key}/=$num;             #frequency of MI
	$cum+=$vmi{$key};
	$vmi{$key}=$cum;              #cumulative frequency of MI
    }
    print $num;
    return %vmi;
    
}
=cut

#parameter: number of bootstrap
sub bootstrap{
    my ($bs_num)=@_;
    my @bpmi;#store the bootstrap probability of each site
    for(my $i=0;$i<$length;$i++){
	for(my $j=0;$j<$length;$j++){
	    $bpmi[$i][$j]=0;
	}   
    }
    my $num=0;#the total number of mi values of all the bootstraps

    
    srand (time ^ $$ ^ unpack "%L*", `ps axww | gzip -f`);
    
    for(my $n=0;$n<$bs_num;$n++){#each bootstrap
    
        my @bseq;
        for(my $j=0;$j<$length;$j++){#each site
            my %rands;
            for(my $i=0;;){
                last if($i==$count);
                my $rn=int(rand($count));#get a seq randomly
                if(exists $rands{$rn}){#if this seq has been used, choose another seq
                    $rands{$rn}=0;
                    next;
                }else{
                    $rands{$rn}=1;
                }
                $bseq[$i++][$j]=$seqs[$rn][$j];#shuffle the bases in one column
                #print $rn." ";
            }
            #print "\n";
        }
    
        my @bmi=&MI(@bseq);#calculate MI of this bootstrap sample
    
	#if in a site, the bootstrap mi value greater than the true mi value, the number add 1 in this site 
        for(my $i=0;$i<$length-1;$i++){
            for(my $j=$i+1;$j<$length;$j++){
                if($bmi[$i][$j]>=$mi[$i][$j]){
			$bpmi[$i][$j]+=1;
                }
            }
        }
        print "$n of the ".$_[0]."\n";
    }
    
    #calculate the distribution of probibality of mi value
    for(my $i=0;$i<$length-1;$i++){
        for(my $j=$i+1;$j<$length;$j++){
	    $bpmi[$i][$j]/=$bs_num;
	    $bpmi[$j][$i]=1;
    	}
	$bpmi[$i][$i]=1;
    }
    $bpmi[$length-1][$length-1]=1;
    
    return @bpmi;
}

