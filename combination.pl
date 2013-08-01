#! perl -w
#generate all the possible combinations of letter 
use strict;

my @ori=('a'..'k');
my @tree=('=');
@tree=addone(@tree);
for my $i (1..11){
	my @top10;
	for my $j (0..9){
		$top10[$j]=$tree[$j];
	}

	for my $t (@top10){
		print $i.$t."\n";
		
	}
	
	@tree=addone(@top10);
	
	
}


sub addone{
	my (@group)=@_;
	my @new=();
	for my $i (@group){
		my @used=split /,/,$i;
		my @left=remove(@used);
		for my $j (@left){
			push @new,$i.",".$j;		
		}
	}
	return unique(@new);
	#return @new;
}

sub remove{
	my (@del)=@_;
	my $end=-1;
	for(my $i=0;$i<=$#ori;$i++){
		if($del[$#del] eq $ori[$i]){
			$end=$i;
			last;
		}
	}
	my @unused=();
	for(my $i=$end+1;$i<=$#ori;$i++){
		push @unused,$ori[$i];
	}
	my @left=();
	for my $i (@unused){
		my $flag=0;
		for my $j (@del){
			if ($i eq $j){
				$flag=1;
				last;
			}
		}
		push @left,$i if($flag==0);
	}
	return @left;
}

sub remove2{
	my (@del)=@_;
	
	my @left=();
	for my $i (@ori){
		my $flag=0;
		for my $j (@del){
			if ($i eq $j){
				$flag=1;
				last;
			}
		}
		push @left,$i if($flag==0);
	}
	return @left;
}

sub unique{
	my(@all)=@_;
	for my $temp (@all){
		$temp=join(",", (sort {$a cmp $b} split /,/,$temp));	
	}
	for(my $i=0;$i<$#all;$i++){
		for(my $j=$i+1;$j<=$#all;$j++){
			if($all[$i] eq $all[$j]){
				$all[$j]="#";
			}
		}
	}
	my @one=();
	for my $temp (@all){
		push @one,$temp if $temp ne "#";
	}
	return @one;
}
