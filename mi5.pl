#! /usr/bin/perl -w
############# weibo 2009-9-21 ####################
############ remove the consistent columns ########
use strict;

use Bio::SeqIO;
my $input=shift;
my $outMI=$input."."."MI3";
my $outMIp=$input."."."MIp";
my $outen=$input."."."En";

my $seqio=Bio::SeqIO->new(-file=>$input,-format=>'fasta');
open(OUT,">$outMI");
open(OUTO,">$outMIp");
open(OUTEN,">$outen");

my $sequence;
my @aa;
while($sequence=$seqio->next_seq){
    my $seq=$sequence->seq;
    push @aa,[split //,$seq];
}

my $row=scalar(@aa);
my $column=scalar(@{$aa[0]});



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

    print OUTEN ($i+1)."\t".$entropy[$t]."\n";
    $t++;
}

$column=scalar(@entropy);


print "The entropy calculation has finished!\n";

print "Now start to calculate MI(a,b)\n";

 my @mi;
 my @MIp;
   
    for(my $i=0;$i<$column;$i++){
	for(my $j=0;$j<$column;$j++){
	    $mi[$i][$j]=0;
	    $MIp[$i][$j]=0;
	}   
    }

my $mmi=0;
    for(my $i=0;$i<$column-1;$i++){
	for(my $j=$i+1;$j<$column;$j++){
	    my %pairs;
	
	    for(my $k=0;$k<$row;$k++){
		my $p=$trim[$k][$i]." ".$trim[$k][$j];
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
        
	    $mi[$i][$j]=($entropy[$i]+$entropy[$j]-$hab);            #...MI(a,b)
            
            $mmi+=$mi[$i][$j];
	}
        
    }
    
$mmi/=$column*($column-1)/2;       # mean MI
print "average MI = $mmi\n";


print "MI calculation has finished!\n";

my @MIax;
for(my $i=0;$i<$column;$i++){
	for(my $j=0;$j<$column;$j++){
	    
	    $mi[$j][$i]=$mi[$i][$j];
	    $MIax[$i]+=$mi[$i][$j];
	   
	    print OUT $mi[$i][$j]."\t";   # output MI array
	}
	$MIax[$i]=($MIax[$i]-$mi[$i][$i])/($column-1);             #...MI(a,x)
        print OUT "\n";
    }

my @apc;


for(my $i=0;$i<$column-1;$i++){
    for(my $j=$i+1;$j<$column;$j++){
        $apc[$i][$j]=$MIax[$i]*$MIax[$j]/$mmi;      #...APC(a,b)
        $MIp[$i][$j]=$mi[$i][$j]-$apc[$i][$j];      #...MIp(a,b)
	$MIp[$i][$j]=0 if$MIp[$i][$j]<0;
        #$MIp[$i][$j]=$apc[$i][$j];
       
    }
    print "...".($i+1)." of $column\n";
}

print "MIp calculation has finished!\n";

for(my $i=0;$i<$column;$i++){
    for(my $j=0;$j<$column;$j++){
        printf OUTO "%2.10f\t",$MIp[$i][$j];                     #print out MIp matrix
        
    }
    print OUTO "\n";

}
