#! perl -w
#实现渐变色
use strict;

open(OUT,">渐变色.html");
print OUT "<html><head><title>渐变色</title><style type=\"text/css\">\n
                .mir{width:100px;height:2px;}</style>\n</head><body>\n";

my @c1=(0,255,255);#start
my @c2=(255,0,0);#end
my $rate;
my($r,$g,$b);
for $r (0..255){
	$rate=($r-$c1[0])/($c2[0]-$c1[0]);
	$g=$rate*($c2[1]-$c1[1])+$c1[1];
	$b=$rate*($c2[2]-$c1[2])+$c1[2];
	#print "($r,$g,$b)\t";
	my $color=sprintf("#%02x%02x%02x",$r,$g,$b);
	print OUT "<div class=\"mir\" style=\"background:".$color."\"></div></a>\n";
}

print OUT "</body></html>";

