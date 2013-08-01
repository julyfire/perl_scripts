#! /usr/bin/perl -w
############# weibo 2009-9-21 ####################
############ remove the consistent columns ########
use strict;

use Bio::SeqIO;
my $input=shift;
my $outMI=$input."."."MI3";
my $outMIp=$input."."."MIp";
my $outbs=$input."."."bs";

my $seqio=Bio::SeqIO->new(-file=>$input,-format=>'fasta');
open(OUT,">$outMI");
open(OUTO,">$outMIp");
open(OUTBS,">$outbs");

my $sequence;
my @aa;
while($sequence=$seqio->next_seq){
    my $seq=$sequence->seq;
    push @aa,[split //,$seq];
}

my $row=scalar(@aa);
my $column=scalar(@{$aa[0]});


my $index;
my @entropy;
my @trim;
my $t=0;

for(my $i=0;$i<$column;$i++){
    my %temp;
    my $p;
    for(my $j=0;$j<$row;$j++){
        $p=$aa[$j][$i];
        if(exists $temp{$p}){
            $temp{$p}+=1;
        }
        else{
            $temp{$p}=1;
        }
    }
    my $ele=keys %temp;
        
    next if($ele==1);
    
    for(my $k=0;$k<$row;$k++){
        $trim[$k][$t]=$aa[$k][$i];
    }
    
    
    while(my($key,$value)=each %temp){
        $entropy[$t]+=-$value/$row*log($value/$row)/log(2);  #...H(a)
            #print "Entropy: ".$key."=>".$value."\n";    
    }
    $t++;
    $index.=$i." ";
    
}

$column=scalar(@entropy);


print "The entropy calculation has finished!\n";

print "Now start to calculate MI(a,b)\n";
my @mi=&MI(@trim);
print "MI calculation has finished!\n";

for(my $i=0;$i<$column;$i++){
	for(my $j=0;$j<$column;$j++){
	    print OUT $mi[$i][$j]."\t";
	}
        print OUT "\n";
    }

print "Now start to bootstrap!\n";
print "Input the bootstrap number: ";
chomp(my $bnum=<STDIN>);

print "Input shuffle number: ";
chomp(my $snum=<STDIN>);

my %fmi=&bootstrap($bnum,$snum);
foreach my $key(sort keys %fmi){
    print OUTBS $key."\t".$fmi{$key}."\n";
}
print "bootstrap calculation has been finished!\n";

sub MI{
    my @seqs=@_;
    my @mi;
   
    for(my $i=0;$i<$column;$i++){
	for(my $j=0;$j<$column;$j++){
	    $mi[$i][$j]=0;
	}   
    }


    for(my $i=0;$i<$column-1;$i++){
	for(my $j=$i+1;$j<$column;$j++){
	    my %pairs;
	
	    for(my $k=0;$k<$row;$k++){
		my $p=$seqs[$k][$i]." ".$seqs[$k][$j];
		if(exists $pairs{$p}){
		    $pairs{$p}+=1;
		}
		else{
		    $pairs{$p}=1;
		}
            
	    }
	    my $hab=0;
	    while(my($key,$value)=each %pairs){                      #...H(a,b)
		$hab+=-$value/$row*log($value/$row)/log(2);
	    }
        
	    $mi[$i][$j]=($entropy[$i]+$entropy[$j]-$hab)/$hab;            #...MI(a,b)/H(a,b)
            
            
	}
        
    }
    
    return @mi;
}

sub bootstrap{
    my %vmi;
    my $num=0;
    my @bseq=@trim;
    
    srand (time ^ $$ ^ unpack "%L*", `ps axww | gzip -f`);
    
    for(my $n=0;$n<$_[0];$n++){
        
        for(my $i=0;$i<$column;$i++){
            for(my $b=0;$b<$_[1];$b++){
                my ($r1,$r2)=(rand($row),rand($row));
                if($r1==$r2){
                    $b-=1;
                    next;
                }
                ($bseq[$r1][$i],$bseq[$r2][$i])=($bseq[$r2][$i],$bseq[$r1][$i]);
                
            }
        
        }
        
        
    
        my @bmi=&MI(@bseq);
    
        for(my $i=0;$i<$column-1;$i++){
            for(my $j=$i+1;$j<$column;$j++){
                my $mii=$bmi[$i][$j];
                if(exists $vmi{$mii}){
                    $vmi{$mii}+=1;
                }else{
                    $vmi{$mii}=1;
                }
                $num++;
            }
        }
        
       print (($n+1)." of the $_[0] bootstrap\n");
    } 
    
    foreach my $key(sort keys %vmi){
        $vmi{$key}/=$num;
    }
    print $num;
    return %vmi;
    
}
