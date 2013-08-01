#! perl -w
#parse the Vienna dot bracket notation format of RNA secondary structure
#remove any basepair the distance of which is greater than a threshold,say 100
#by weibo 2010-7-12
use strict;

my($seq,$str)=('','');
if(@ARGV){
  my $infile=shift;
  open(IN,$infile) or die "cannot open the file\n";
  $seq=<IN>;
  $str=<IN>;
  $str=~s/ \(.*\)$//;
}
if($str eq ''){
  print "Input a sequence string:\n>";
  $seq=<STDIN>;
  print "Input a dot bracket string:\n>";
  $str=<STDIN>;
}
chomp($seq);
chomp($str);
my @loc;
my @base=split //,$str;

print "set a threshold:\n>";
my $threshold=<STDIN>;
chomp($threshold);

my($i,$j,$c,$value,$sum,$lo)=(0,0,0,0,0,0);
for($i=0;$i<$#base;$i++){
  next unless $base[$i] eq '('; #find first left brackets
  for($j=$i+1;$j<=$#base;$j++){
    next unless $base[$j] eq ')'; #find first right brackets
    $sum=0;
    for($c=$i;$c<=$j;$c++){
      $value=0;
      $value=1 if $base[$c] eq '(';
      $value=0 if $base[$c] eq '.';
      $value=-1 if $base[$c] eq ')';
      $sum+=$value;
    }
    next unless $sum==0;
    if($j-$i>=$threshold){
      $base[$i]='.';
      $base[$j]='.';
    }
    else{
      $loc[$lo++]=[$i,$j];
      $i=$j;
    }
    last; # end the inner loop when finding a local structure
  }
}
$str=join '',@base;
print "result:\n".$str."\n";

my $e=0;
for my $f (@loc){
#print "start:".$f->[0].", end:".$f->[1]."\n";
my $seg=substr($seq,$$f[0],$$f[1]-$$f[0]+1);
my $segstr=substr($str,$$f[0],$$f[1]-$$f[0]+1);
print $seg."\n";
print $segstr."\n";
$e+=`./dG '$seg' '$segstr'`;
}
print "energy=".$e."\n";
#-------------------------------------------------------------


