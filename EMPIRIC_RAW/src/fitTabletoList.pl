#!/usr/bin/perl


$inputTable = $ARGV[0];
#($name, $ext) = split (/\./, $inputTable);
($name, $ext) = $inputTable =~ m/(.*)\.(.*)$/;

print "-- $name -- $ext --\n";

$outputList = $name."-List.".$ext;

print "$outputList\n";
#exit(1);

open IN, $inputTable or die "cannot open input\n"; 

$line = <IN>; #remove first line 

while(<IN>)
{
	chomp;
	my @mut = split (/,/, $_);
	shift(@mut);
	push @fitness,[@mut];  
}
close IN;

open OUT, "> $outputList" or die "cannot open output file\n";

for $i ( 0 .. $#fitness ) {
        for $j ( 0 .. $#{$fitness[$i]} ) {
		$f = $fitness[$i][$j]; 
		if ($f != 0){
		$newF = sprintf("%.3f", $f);
            	print "$j $i $newF\n";
	    	print OUT "$j $i $newF\n";	
		}
        }
}
close OUT; 
