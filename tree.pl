#! perl -w
use strict;
use Bio::TreeIO;
my $infile=shift;
my $input=new Bio::TreeIO(-file=>$infile,
			-format=>"newick",
			-internal_node_id=>"bootstrap");
while(my $tree=$input->next_tree){
	$tree->move_id_to_bootstrap;
	my %dis;
	my $max=0;
	my $refseq=$tree->find_node(-id=>"GU131785");
	for my $taxa($tree->get_leaf_nodes){
		my $d=$tree->distance(-nodes=>[$refseq,$taxa]);
		$dis{$taxa->id}=$d;
		if($max<$d){
			$max=$d;
		}
	}
	
	#foreach my $key(sort {$dis{$b}<=>$dis{$a}} keys %dis){
	#	print $key.":".$dis{$key}."\n";
	#}
	
	my $interval=$max/100;
	my @group;	
	for(my $i=0;$i<=$max-$interval;$i+=$interval){
		my %t;
		foreach my $key(sort {$dis{$b}<=>$dis{$a}} keys %dis){
			if($i<$dis{$key} && $dis{$key}<=$i+$interval){
				$t{$key}=$dis{$key};
			}			
		}
		push @group,\%t;
	}
	
	foreach my $item(@group){
		next unless %$item;
		my $key=shift @{[sort {$dis{$b}<=>$dis{$a}} keys %$item]};
		print $key."\t".$item->{$key}."\n";
		#print $key."\n";
	}
	
	
}

