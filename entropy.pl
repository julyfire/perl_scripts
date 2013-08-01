#! /usr/bin/perl -w

############ weibo 2009-9-9 14:39 #########################
############ mi & apc (original version, so slow!) #########

use strict;
use Bio::SeqIO;
my $input=shift;
my $outMI=$input."."."MI";
my $outMIp=$input."."."MIp";
my $temp=$input."."."temp";

my $seqio=Bio::SeqIO->new(-file=>$input,-format=>'fasta');
open(OUT,">$outMI");
open(OUTO,">$outMIp");
open(TEMP,">$temp");

my $length=$seqio->next_seq->length;

my $key;
my $value;


####################################     figure out the entropy
my @entropy;

print TEMP "#################### entropy ####################\n";
for(my $i=0;$i<$length;$i++){
    my %sites;
    
    my $sequence;
    my $count=0;
    $seqio=Bio::SeqIO->new(-file=>$input,-format=>'fasta');
    while($sequence=$seqio->next_seq){
        my $seq=$sequence->seq;
        my $p=substr($seq,$i,1);
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
    print "...".($i+1)." of $length\n";
    print TEMP "site ".$i."'s entropy is ".$entropy[$i]."\n";
}
print "The entropy calculation has finished!\n";

########################################     figure out the MI
print "Now start to calculate MI(a,b)\n";
my @mi;
my @MIp;
my @apc;
for(my $i=0;$i<$length;$i++){
    for(my $j=0;$j<$length;$j++){
        $mi[$i][$j]=0;
        $MIp[$i][$j]=0;
        $apc[$i][$j]=0;
    }
}

my $MIm=0;

print TEMP "\n############# H(a,b) ###############\n";
for(my $i=0;$i<$length-1;$i++){
        
    for(my $j=$i+1;$j<$length;$j++){
        my %pairs;
        $seqio=Bio::SeqIO->new(-file=>$input,-format=>'fasta');
        my $count=0;
        while(my $sequence=$seqio->next_seq){
            my $seq=$sequence->seq;
            my $site1=substr($seq,$i,1);
            my $site2=substr($seq,$j,1);
            my $p=$site1."-".$site2;
            if(exists $pairs{$p}){
                $pairs{$p}+=1;
            }
            else{
                $pairs{$p}=1;
            }
            $count++;
        }
        my $hab=0;
        while(($key,$value)=each %pairs){                      #...H(a,b)
            $hab+=-$value/$count*log($value/$count)/log(20);
        }
        
        print TEMP "Hab: pair($i,$j)=".$hab."\n";
        $mi[$i][$j]=$entropy[$i]+$entropy[$j]-$hab;            #...MI(a,b)
       
        #print TEMP "MIab: pair($i,$j)=".$mi[$i][$j]."\n";
        $MIm+=$mi[$i][$j];
    }
    print "...".($i+1)." of $length\n";
  
}

print "MI calculation has finished!\n";


##########################################     figure out the MIp
print "Now start to calculate MIp(a,b):\n";

print TEMP "\n############# mean MI ############3\n";
$MIm=$MIm*2/($length*($length-1));                             #...mean MI
print TEMP "mean MI=$MIm\n";

my @MIax;

print TEMP "\n################ MI(a,x) ###############\n";
for(my $i=0;$i<$length;$i++){
    for(my $j=0;$j<$length;$j++){
        printf OUT "%2.10f\t",$mi[$i][$j];                     #print out MI matrix
        $mi[$j][$i]=$mi[$i][$j];
        $MIax[$i]+=$mi[$i][$j];
        
    }
    print OUT "\n";
    $MIax[$i]=($MIax[$i]-$mi[$i][$i])/($length-1);             #...MI(a,x)
    print TEMP "MI($i,x)=$MIax[$i]\n";
}

print TEMP "\n#################### MI/APC #########################\n";

for(my $i=0;$i<$length-1;$i++){
    for(my $j=$i+1;$j<$length;$j++){
        $apc[$i][$j]=$MIax[$i]*$MIax[$j]/$MIm;      #...APC(a,b)
        $MIp[$i][$j]=$mi[$i][$j]-$apc[$i][$j];      #...MIp(a,b)
        
        print TEMP "MI/APC($i,$j)=$mi[$i][$j]/$apc[$i][$j]=".$mi[$i][$j]/$apc[$i][$j]."\n" if($apc[$i][$j]!=0);
    }
    print "...".($i+1)." of $length\n";
}

print "MIp calculation has finished!\n";

for(my $i=0;$i<$length;$i++){
    for(my $j=0;$j<$length;$j++){
        printf OUTO "%2.10f\t",$MIp[$i][$j];                     #print out MIp matrix
        
    }
    print OUTO "\n";

}





