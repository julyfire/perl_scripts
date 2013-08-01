#! /usr/bin/perl -w
use strict;
my @aa=('a'..'z');
 srand (time ^ $$ ^ unpack "%L*", `ps axww | gzip -f`);

my %rands;
for(my $i=0;;){
    last if($i==26);
    my $rn=int(rand(26));
    next if(exists $rands{$rn});  
    print $aa[$rn]."\n";
    $i++;
}

