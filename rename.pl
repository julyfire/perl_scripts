#! perl -w
#replace the special character of filenames with "_"
#2010-6-9 wb 
use strict;
#method 1:
#my @all=glob ".* *";
#foreach my $fn (@ARGV){}

my $dir=shift;
$dir=~s/\/$//;
#method 2:
#my @all=<$dir/*>;
#foreach my $fn (@all){}
#method 3:
opendir(DH,$dir);
foreach my $fn (readdir DH){
   next unless($fn=~/[\s\n\r\:\|\!\?\"]+/);
   (my $newfn=$fn)=~s/\s|\n|\?/_/g;
   if(-e $newfn){
     warn "can't rename $fn to $newfn: $newfn exists\n";
   }
   elsif(rename $dir."/".$fn, $dir."/".$newfn){
   }else{
      warn "rename $fn to $newfn failed: $!\n";
   }
}
closedir DH;

