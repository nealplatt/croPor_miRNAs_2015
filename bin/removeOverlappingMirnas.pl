#takes in merged bed file and finds best scoreing

#bedtools merge -d -1 -i BEDFILE | bedtools intersect -wao -a - -b BEDIFILE | perl removeOverlappingMirnas.pl


while(<>){
	chomp;
	($M_scaffold, 	$M_start, 	$M_end, 	$M_name, 	
	$M_e, 		$M_f, 		$name,		$score, 
	$orient,	$start, 	$stop, 		$color, 	
	$overlap) = split ("\t", $_);

	

		$uniqID="$M_scaffold.$M_start.$M_end";
		
		$list{$uniqID}++;

		if( ($list{$uniqID} > 1) && ($score <= $score{$uniqID}) ){
				
			next;
		}

		$name{$uniqID} = $name;		
		$score{$uniqID} = $score; 
		$start{$uniqID} = $start;	
		$stop{$uniqID} = $stop;
		$orient{$uniqID} = $orient;
		$scaffold{$uniqID} = $M_scaffold;
}



foreach $key (keys (%list)){
		
	
		print "$scaffold{$key}\t$start{$key}\t$stop{$key}\t$name{$key}\t$score{$key}\t$orient{$key}\n";
}
