#!/usr/bin/perl


@refseqnames = split(",", $ARGV[0]); # b3,b4 comma-separated list of library names

$basedir = "../output";
$outdir = "../output";

$finaloutput = $ARGV[1];

foreach $r (@refseqnames) {

$fnlist .= "$basedir/fit-table-$r-aa.dat ";

}

$fnout = "$outdir/concat-table.temp";
$cmd =  "paste $fnlist > $fnout";
system($cmd);

open FOUT, "> $finaloutput";

open F, "$fnout";
while(<F>) {
$line =  $_;

$line =~ s/\t/,/g;
$line =~ s/\s//g;

#print $line;

@cols = split(",", $line);

$rewrittenline = "";
$i = 0;
foreach $c (@cols) {
	#print "$i ";
	if (!($i%11 == 0 && $i!=0)) { $rewrittenline .= "$c," };
	#print "\n";
	$i++;
	}

$rewrittenline =~ s/,$//;


print FOUT "$rewrittenline\n";

}


close FOUT;

