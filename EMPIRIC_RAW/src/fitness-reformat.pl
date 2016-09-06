#!/usr/bin/perl

$fitnessfn = $ARGV[0];
$outputfn = $ARGV[1];
$betaout = $ARGV[2];

$numpositions = 10;
$numelements = 21; # A..Y*

# input is ASSUMED to be in this order:
$aaorder_old = "ACDEFGHIKLMNPQRSTVWY*";

# output order: 
$aaorder_new = "ILVFYWCMAGNSTHQDEKRP*";


open F, "$fitnessfn";
for($j=0; $j<$numpositions; $j++){

	for($i=0; $i<$numelements;$i++) {

		$line[$i] = <F>;
		@f = split(" " , $line[$i]);

#		if ($f[2] < -1.5) { $f[2]=-1.5; }
#		if ($f[2] > 0.5 ) { $f[2]= 0.5; }


		$fitness[$i][$j] = $f[2];
		if ($i == 0) {
			($pos[$j]) = $f[1] =~ m/(\D\d+)\D/;
			($posnumber[$j]) = $f[1] =~ m/\D(\d+)\D/;
			}

		($letter[$i]) = $f[1] =~ m/.*(\D)/;

		}

}
close F;


for($j=0; $j<$numpositions; $j++){

for($i=0; $i<$numelements; $i++) {

$aaname = substr($aaorder_old, $i, 1);

for($k=0;$k<$numelements;$k++) { if (substr($aaorder_new, $k,1) eq $aaname) { $newindex = $k; last; } } 

$fitness_new[$k][$j]=$fitness[$i][$j];

}


}


open OUT, "> $outputfn";
print OUT "XXX,";
for($j=0;$j<$numpositions;$j++) { 
	print OUT "$pos[$j]"; 
	if ($j!=$numpositions-1) { print OUT ","; }
	}
print OUT "\n";

for($i=0;$i<$numelements;$i++) {
	
	print OUT substr($aaorder_new,$i,1).",";

	for($j=0;$j<$numpositions;$j++) {
		print OUT "$fitness_new[$i][$j]";
		if ($j!=$numpositions-1) { print OUT ","; }
		}
	print OUT "\n";
	}
close OUT;

# beta factors for the PDB
open OUT, "> $betaout";
for($j=0; $j<$numpositions; $j++){
$s = 0;
for($i=0; $i<$numelements-1; $i++) {
$s += $fitness[$i][$j];
}
$s /= ($numelements-1);
print OUT "$posnumber[$j] $s\n";
}
close OUT;

