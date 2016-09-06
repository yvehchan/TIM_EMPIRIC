#!/usr/bin/perl
use Switch;

@refseqnames = split(",", $ARGV[0]); # b3,b4 comma-separated list of library names

$referencefastafn = $ARGV[1];

@codonlist =  ( "AAA","AAT","AAG","AAC",
"ATA","ATT","ATG","ATC",
"AGA","AGT","AGG","AGC",
"ACA","ACT","ACG","ACC",
"TAA","TAT","TAG","TAC",
"TTA","TTT","TTG","TTC",
"TGA","TGT","TGG","TGC",
"TCA","TCT","TCG","TCC",
"GAA","GAT","GAG","GAC",
"GTA","GTT","GTG","GTC",
"GGA","GGT","GGG","GGC",
"GCA","GCT","GCG","GCC",
"CAA","CAT","CAG","CAC",
"CTA","CTT","CTG","CTC",
"CGA","CGT","CGG","CGC",
"CCA","CCT","CCG","CCC" );

@mmei_sites = ("TCCAAC", "TCCGAC",   "GTTGGA", "GTCGGA" );


foreach $refseqname (@refseqnames) {

print "$refseqname\n";

$outputdir = "../output";
#$inputdir = "../input";

$refseqfn = "$referencefastafn";
$mapfn =  "$outputdir/maptable-$refseqname.txt";

$outputfn = "$outputdir/mmeisites-$refseqname.txt";

open OUTF, "> $outputfn";

open F, $refseqfn;
while(<F>) {
chomp;
($f1,$f2) = split(" ", $_);
if ($f1 eq $refseqname) { $refseq = $f2; }
}
close F;

@mmeisites=();
@mmeicodons=();

#b3 0 S TCT 104 S
open F, "$mapfn" || die "cannot open $mapfn\n";
while(<F>) {
chomp;
@f = split(" ", $_);
$pdb_resnum[$f[1]] = $f[4];
}
close F;




for($i=0;$i<10;$i++) {

foreach $mutantcodon (@codonlist) {

$refseqmutant = $refseq;

substr($refseqmutant, $i*3, 3) = $mutantcodon;

#print "$refseqmutant\n";

# look for MmeI site in the sequence:
$flag = 0;
foreach $site (@mmei_sites) {
if ( index($refseqmutant,$site)>-1 ) { $flag = 1;}
}

if ($flag == 1) { #print "codon $mutantcodon at pos $i creates an MmeI site\n"; 
		push @mmeisites, $i;
		push @mmeicodons, $mutantcodon;
		} 

}

}

for($i=0;$i<=$#mmeisites;$i++) {

$wtcodon = substr($refseq, $mmeisites[$i]*3,3);

$wt_aa  = &getAA($wtcodon);
$mut_aa = &getAA($mmeicodons[$i]);

print "codon $mmeicodons[$i] at pos $mmeisites[$i] creates an MmeI site, $wt_aa $pdb_resnum[$mmeisites[$i]] $mut_aa\n";

print OUTF "$refseqname $wt_aa$pdb_resnum[$mmeisites[$i]]$mut_aa\n";

}

close OUTF;

}



sub getAA($AA){ 
my $codon = @_[0];
switch ($codon){
case [ATT, ATC, ATA] {$aa = "I";}
case [CTT, CTC, CTA, CTG, TTA, TTG] {$aa = "L";}
case [GTT, GTC, GTA, GTG] {$aa = "V";}
case [TTT, TTC] {$aa = "F";}
case [ATG] {$aa = "M";}
case [TGT, TGC] {$aa = "C";}
case [GCT, GCC, GCA, GCG] {$aa = "A";}
case [GGT, GGC, GGA, GGG] {$aa = "G";}
case [CCT, CCC, CCA, CCG] {$aa = "P";}
case [ACT, ACC, ACA, ACG] {$aa = "T";}
case [TCT, TCC, TCA, TCG, AGT, AGC] {$aa = "S";}
case [TAT, TAC] {$aa = "Y";}  
case [TGG] {$aa = "W";}
case [CAA, CAG] {$aa = "Q";}
case [AAT, AAC] {$aa = "N";}
case [CAT, CAC] {$aa = "H";}
case [GAA, GAG] {$aa = "E";}
case [GAT, GAC] {$aa = "D";}
case [AAA, AAG] {$aa = "K";}
case [CGT, CGC, CGA, CGG, AGA, AGG] {$aa = "R";}
case [TAA, TAG, TGA] {$aa = "*";} 
}
return $aa; 
}




