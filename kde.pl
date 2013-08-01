#! perl -w
use strict;
use Statistics::KernelEstimation;


my $infile=shift;
open(IN,$infile);
my @d=();
my @label=();
my $i=0;
my $j=0;
while(<IN>){
	chomp;
	s/,,/,/;
	my @line=split /,/;
	my $name=shift @line;
	push @label,$name;
	for($j=0;$j<scalar(@line);$j++){
		$d[$i][$j]=$line[$j];
		$d[$j][$i]=$line[$j];
	}
	$i++;
}
close IN;

my @data=();
for($i=0;$i<scalar(@label)-1;$i++){
	for($j=$i+1;$j<scalar(@label);$j++){
		push @data,$d[$i][$j];
	}
}

my $s = Statistics::KernelEstimation->new();

for my $x ( @data ) {
    $s->add_data( $x );
}

my $w = $s->default_bandwidth();
my ( $min, $max ) = $s->extended_range();

print("default:".$w."\n"."fast:".$s->optimal_bandwidth()."\n"."safe:".$s->optimal_bandwidth_safe()."\n");

my %result;
for( my $x=$min; $x<=$max; $x+=($max-$min)/100 ) {
	$result{$x}=$s->pdf( $x, $w );
    	print $x, "\t", $s->pdf( $x, $w ), "\n";
}
my @keys=sort{$a<=>$b} keys %result;
for(my $i=1;$i<$#keys;$i++){
	if($result{$keys[$i]}<$result{$keys[$i-1]} and $result{$keys[$i]}<$result{$keys[$i+1]}){
		print $keys[$i]."\n";
		last;
	}
}
	

