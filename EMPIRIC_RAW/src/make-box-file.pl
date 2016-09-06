#!/usr/bin/perl

# creates box files for the heatmap

$inputtable = $ARGV[0];  # eg  fit-table-b3b4-aa.dat

@warnfiles = split(",", $ARGV[1]); # b3,b4 comma-separated list of library names

$outputdir = "../output";

# black boxes around WT

open F, $inputtable;
$firstline = <F>;
chomp($firstline);
while(<F>) {
	@f = split(",", $_);
	$aatype[$rowcount++] = $f[0];
	}
close F;

@wtcode = split(",",$firstline);
for($i=1;$i<=$#wtcode;$i++) {
	($wtaa) = $wtcode[$i] =~ m/^(\D).*/;
	($position[$i]) = $wtcode[$i] =~ m/^\D(\d+)/;
	for($j=0;$j<$rowcount;$j++) { if ($wtaa eq $aatype[$j]) { $row = $j; } }
	$col = $i-1;
	print "$col $row k square\n";
	}

# sites with zero abundance

foreach $fn (@warnfiles) {

open F, "$outputdir/slopes-warn-$fn-aa.dat";
while(<F>) {
chomp;
($lib, $wt, $pos, $mut ) = $_ =~ m/(.*)\s+(\D)(\d+)(\D)/;

for($j=1;$j<$#wtcode;$j++) { if ($pos eq $position[$j]) { $col = $j-1; } }
for($j=0;$j<$rowcount;$j++) { if ($mut eq $aatype[$j]) { $row = $j; } }
	print "$col $row r squarecirc\n";
}
close F;

}

# MmeI sites

foreach $fn (@warnfiles) {

open F, "$outputdir/mmeisites-$fn.txt";
while(<F>) {
chomp;
($lib, $wt, $pos, $mut ) = $_ =~ m/(.*)\s+(\D)(\d+)(\D)/;

for($j=1;$j<$#wtcode;$j++) { if ($pos eq $position[$j]) { $col = $j-1; } }
for($j=0;$j<$rowcount;$j++) { if ($mut eq $aatype[$j]) { $row = $j; } }
	print "$col $row r squarecross\n";
}
close F;

}




# TODO: sites with too low abundance in the initial library
