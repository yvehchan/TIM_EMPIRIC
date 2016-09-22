#!/share/bin/perl
use Switch; 
use Text::Tabs;
use Term::ANSIColor; 

### This script will generate primers and cassettes necessary to generate the EMPIRIC libraries
### Requires the DNA sequence, first residue of fasta sequence and first residue to begin 10 position randomization  
###example of usage### 
##./DeepSeqPrimers.pl SsIGPS-noTag-AA27-248.fa 27 104

$tabstop = 4; 

$BSAI = "ggtctc"; #cut site for cloning in cassette
$SPHI = "gcatgc"; #cut site to prevent religation 
$filler = "a"; #filler to remove offset of BSA restriction site
$five_primer = "ggtggtggt"; #random DNA seq, will be cut out 
$three_primer = "gcagcagca"; #randowm DNA seq, will be cut out

my $seqFile = $ARGV[0]; #fasta sequence, remove tags
my $gene = substr $seqFile, 0, 2; #first five letters of file name
my $seqStart = $ARGV[1]; #first residue of fasta sequence 
my $AAstart = $ARGV[2]; #residue to begin 10 position randomization 
$AAend = $AAstart+9; #inclusive of ends 
$start = $AAstart-7; #display preceding sequence
$end = $AAend+7;  #display succeeding sequence 
$ForStart = $AAend+2; #end cassette 2 aa after randomization 
$ForEnd = $ForStart+7; #use 7 aa length for primer annealing 
$RevStart = $AAstart-3; #start cassette 2 aa before randomization
$RevEnd = $RevStart-7; #use 7 aa length for primer annealing
@ForPrimer = ();
$fPrimer = "";
@RevPrimer = ();
$rPrimer = "";

#print "Start: ".$start."\t\tEnd: ".$end."\n"; 
print "Randomization start: ".$AAstart."\tRandomization end: ".$AAend."\n";
print "Cassette start: ".($RevStart+1)."\t\tCassette end: ".$ForStart."\n";

open IN, $seqFile or die;
$DNAseq = <IN>;
close IN;

$order = "$gene$AAstart\-$AAend"."EMPIRIC_Primers.txt";
$order2 = "$gene$AAstart\-$AAend"."Illumina_Primers.txt";

open (OUT, "> $order"); 
open (OUT2, ">$order2");

@DNA3 = ($DNAseq =~ m/.../g); #every 3 DNA position sequence 
for ($i = 0; $i < $seqStart-1; $i++){
unshift(@DNA3, 0); #corrects index of fasta sequence from first residue read
}
$size =  scalar (@DNA3); #Number of aa's in sequence

for ($i = $start; $i <=$end; $i++){
       print expand ("$i\t");
}print "\n"; 
for ($i = $start; $i <=$end; $i++){
	$AA = &getAA($DNA3[$i-1]); 
        print expand ("$AA\t");
}print "\n"; 
for ($i = $start; $i <=$end; $i++){
        print expand ("$DNA3[$i-1]\t");
}print "\n";
for ($i = $start; $i <=$end; $i++){
	$comp = &getComp($DNA3[$i-1]);
	print expand ("$comp\t"); 
}print "\n";

for ($i = $ForStart; $i <$ForEnd; $i++){
	my $string = $DNA3[$i];
	push (@ForPrimer, $string); 
}
$fseq = join("", @ForPrimer);
$tempfoverhang3 = substr($fseq, 0, 4);
$foverhang3 = &getComp($tempfoverhang3); 
$fPrimer = join "", $five_primer, $SPHI, $BSAI, $filler, $fseq;
$fPrimerComp = &getComp($fPrimer);

#print $fPrimerComp."\n"; 
print "\nRestriction sequences: SPHI - ";
#print colored ("$SPHI", 'red'); print  "\n";
print "$SPHI\n";
print "Restriction sequences: BSAI - ";
#print colored ("$BSAI", 'yellow'); print "\n";
print "$BSAI\n";
print "\nWhole Primer Forward:\t"; 
print $fPrimer."\n";
print OUT "$gene$AAstart\-$AAend\_VecF\t$fPrimer\n"; 
#print $five_primer;  print colored ("$SPHI", 'red'); print colored ("$BSAI", 'yellow'); print colored ("$fseq", 'green'); print "\n";
print "Whole Primer F-Comp:\t"; print $fPrimerComp."\n"; print "\n"; 
 
my @temp = ();
for ($i = $RevStart-1; $i <=$RevStart; $i++){
        my $string = $DNA3[$i-1]; 
#print "$string"; 
        push (@temp, $string);
}
my $tempOverHang = join("", @temp); 
#print "---------$tempOverhang";
$rOverhang5 = substr($tempOverHang, 2);

#print "****$rseq\n"; 
#rseq = join("", @RevPrimer);

for ($i = $RevStart; $i >$RevEnd; $i--){
        my $string = $DNA3[$i-1];
        $revCompSeq = &getRevComp($string);
        push (@RevPrimer, $revCompSeq);
}
$rseq = join("", @RevPrimer);

#print "----$rseq---\n"; 
$rPrimer = join "", $three_primer,$SPHI,$BSAI, $filler, $rseq; 
$rPrimerComp = &getComp($rPrimer);
print "Whole Primer Reverse:\t"; 
#print $three_primer; print colored ("$SPHI", 'red'); print colored ("$BSAI", 'yellow'); print colored ("$rseq", 'green'); print "\n";
print "$rPrimer\n"; 
print "Whole Primer R-Comp:\t"; print $rPrimerComp."\n";
print OUT "$gene$AAstart\-$AAend\_VecR\t$rPrimer\n";

print "\nCassette overhangs: $rOverhang5----\n";
print "                        ----$foverhang3\n";
print "\n14 aa cassette randomized at 10 positions\n"; 

#----Cassette----
@doubles =();
my $revOverhang = &getRev($foverhang3);

$preRand = $DNA3[($AAstart-3)].$DNA3[($AAstart-2)]; #2 aa before randomization
$postRand = $DNA3[($AAend)].$DNA3[($AAend+1)]; #2 aa after randomization 
for ($i = ($AAstart); $i <=($AAend); $i++){
@cassette = ();
for ($j = $AAstart; $j <=$AAend; $j++){
	$string = $DNA3[$j-1];
 	if ($j == $i){
		$tempPrint ="";
		$AA = &getAA($DNA3[$i-1]);
		print "Randomize: $AA$i\n";
		$tempPrint = "$gene\_"."$AA"."$i";
		$string = "NNN"; 	
	}
        push (@cassette, $string);
	}
	my $casSeqTemp = join("", @cassette);
	my $casSeq = $preRand.$casSeqTemp.$postRand; 
	my $revCasSeq = &getRevComp($casSeq);
	#my $revOverhang = &getRev($foverhang3);
	push (@doubles, $casSeq); 
	print OUT "$tempPrint"."_F"."\t$rOverhang5$casSeq\n";
	print OUT "$tempPrint"."_R"."\t$revOverhang$revCasSeq\n";
	print "Forward: "; 
	print "$rOverhang5$casSeq\n";  
	print "Reverse: ";
	print "$revOverhang$revCasSeq\n";
        print "\n";
}

$sizeDoubles = scalar(@doubles);
for ($i=0; $i<$sizeDoubles-2; $i++){
	$j = $i+2;
	$doubleSeq = $doubles[$i];
	$position = $j*3; 
	substr($doubleSeq, $position, 3) = "NNN"; 
	$revDouble = &getRevComp($doubleSeq);
	$AA = &getAA($DNA3[$i+$AAstart-1]);
	$AA2 = &getAA($DNA3[$j+$AAstart-1]);
	$pos1 = $i+$AAstart;
	$pos2 = $j+$AAstart; 
	print "Randomize: $AA$pos1 and $AA2$pos2\n";
        print "Forward: ";
	print "$rOverhang5$doubleSeq\n";
        print "Reverse: ";
	print "$revOverhang$revDouble\n\n";
}

#------PCR and Adaptor ---
$PCRfiller = "gggaccaccacc";
$MMEI = "tccgac"; 
$Illumina5prime = "acacgacgctcttccgatct";
$Illumina3prime = lc("GGCATACGAGCTCTTCCGATC");

my @temp3 = ();
for ($i = $AAstart-7; $i <$AAstart; $i++){
        my $string = $DNA3[$i-1];
        push (@temp3, $string);
}

$PCR5recog = join("", @temp3);
$PCR5recog = substr($PCR5recog, 0, 20);

my @temp4 = (); 
for ($i = $AAend+48; $i > $AAend+41; $i--){ #123 nt or 41 aa  away from last randomized region
        my $string = $DNA3[$i-1];
	my $string2 = &getRevComp($string); 
        push (@temp4, $string2);
}

$PCR3recog = join("", @temp4);

print "PCR1 to add 5' MMEI and 3' Illumina universal primer sequence to the 3' end\n"; 
print "F: 5' $PCRfiller$MMEI$PCR5recog 3'\n"; #MME1 site cuts 20 nt away on 5' recognition site and 18 nt on 3', leaves 2 nt overhang
print "R: 5' $Illumina3prime$PCR3recog 3'\n\n";

print OUT2 "PCR1 to add 5' MMEI and 3' Illumina universal primer sequence to the 3' end\n";
print OUT2 "F: 5' $PCRfiller$MMEI$PCR5recog 3'\n"; #MME1 site cuts 20 nt away on 5' recognition site and 18 nt on 3', leaves 2 nt overhang
print OUT2 "R: 5' $Illumina3prime$PCR3recog 3'\n\n";


$AdaptorOverhang = substr($DNA3[$AAstart-2], 0, 2);
@Barcodes = (AAA, ACC, AGG, ATT, CAC, CCG, CGT, CTA, GAG, GCT, GGA, GTC, TAT, TCA, TGC, TTG); 
$template = $Illumina5prime."NNN".$AdaptorOverhang;
print "Adaptors:\n";
print OUT2 "Adaptors:\n";
foreach (@Barcodes){
	substr ($template, length($Illumina5prime), length($_), $_);
	print "Barcode: $_ \n";
	print "5' Adaptor: 5' ".$template." 3'\n"; #prints out 5' F adaptors with barcodes
	print "3' Adaptor: 5' ".&getRevComp($_).&getRevComp($Illumina5prime)." 3'\n\n";

	print OUT2 "Barcode: $_ \n";
        print OUT2 "5' Adaptor: 5' ".$template." 3'\n"; #prints out 5' F adaptors with barcodes
        print OUT2 "3' Adaptor: 5' ".&getRevComp($_).&getRevComp($Illumina5prime)." 3'\n\n";

}
	  

close OUT;
close OUT2; 

#-----
sub getRev($string){

my $codon=""; my $rev="";
$codon=@_[0];
#print "getRev $codon\t"; 
#print $codon."\t"; 
$rev = reverse $codon;
return $rev; 
}

#-----
sub getRevComp($string){
my $codon=""; 
my $rev="";
my $revComp="";

$codon=@_[0];
#print "---$codon\t"; 
$rev = &getRev($codon); 
#print $rev."\t"; 
$revComp = &getComp($rev); 
#print $revComp."\n"; 
return $revComp; 
}

#-----
sub getComp($string){
my $codon=@_[0]; 
#print $codon."\n"; 
my @array =();
my @comp =();

@array = split ("", $codon);

foreach(@array){
switch($_){
case [A]{$comp = "T"}
case [T]{$comp = "A"}
case [C]{$comp = "G"}
case [G]{$comp = "C"}
case [a]{$comp = "t"}
case [t]{$comp = "a"}
case [c]{$comp = "g"}
case [g]{$comp = "c"}
case [N]{$comp = "N"}
}
push (@comp, $comp);  
}
$compCodon = join("", @comp); 
return $compCodon; 
}

#-----
sub getAA($AA){
my $codon = @_[0]; 
switch ($codon){
case [ATT, ATC, ATA] {$aa = "I"}
case [CTT, CTC, CTA, CTG, TTA, TTG] {$aa = "L"}
case [GTT, GTC, GTA, GTG] {$aa = "V"}
case [TTT, TTC] {$aa = "F"}
case [ATG] {$aa = "M"}
case [TGT, TGC] {$aa = "C"}
case [GCT, GCC, GCA, GCG] {$aa = "A"}
case [GGT, GGC, GGA, GGG] {$aa = "G"}
case [CCT, CCC, CCA, CCG] {$aa = "P"}
case [ACT, ACC, ACA, ACG] {$aa = "T"}
case [TCT, TCC, TCA, TCG, AGT, AGC] {$aa = "S"}
case [TAT, TAC] {$aa = "Y"}
case [TGG] {$aa = "W"}
case [CAA, CAG] {$aa = "Q"}
case [AAT, AAC] {$aa = "N"}
case [CAT, CAC] {$aa = "H"}
case [GAA, GAG] {$aa = "E"}
case [GAT, GAC] {$aa = "D"}
case [AAA, AAG] {$aa = "K"}
case [CGT, CGC, CGA, CGG, AGA, AGG] {$aa = "R"}
case [TAA, TAG, TGA] {$aa = "*"} 
}
return $aa; 
}
