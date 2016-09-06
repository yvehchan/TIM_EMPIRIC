#!/usr/bin/perl

use Switch;

$tablename     = $ARGV[0]; # outtable-XX-by-time.txt
$numtimepointsfromcmdline = $ARGV[1];
$datatype      = $ARGV[2]; # nt | aa | codons
$sample_to_plot = $ARGV[3];
$slopefn       = $ARGV[4]; # slope file for plotting fit data
$mapfn         = $ARGV[5];  # position mapping to PDB file
$ctrlsample    = $ARGV[6];  # control sample (b3-lib etc) for noise levels
$blacklistfn   = $ARGV[7];  # blacklisted samples/timepoints


$plot_fitted_lines = 1;

$tempdir = "../temp";
$outputdir = "../output";
$plotdir = "../plots";

$outfnbase = "abundance-$sample_to_plot-$datatype-grouped";
$freqfn = "$outputdir/freq-$ctrlsample-$datatype.dat";


$datatype = "codons"; 


$aalist = "ACDEFGHIKLMNPQRSTVWY*";

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


#if ( $datatype eq "nt" )     { $numelements = 4;}
#if ( $datatype eq "aa" )     { $numelements = 21;}
if ( $datatype eq "codons" ) { $numelements = 64;}

# get doubling time for rescaling x axis tics
$dbltimefn = "$outputdir/dbltime-$sample_to_plot.dat";
open F, "$dbltimefn";
$doublingtime = <F>;
close F;

$blacklistlen=0;
open F, "$blacklistfn";
while(<F>) {
chomp;
push @blacklist, $_;
}
close F;  
foreach $x (@blacklist) { print "BLACKLISTED: $x\n"; }


#b3 0 S TCT 104 S
open F, "$mapfn" || die "cannot open $mapfn\n";
while(<F>) {
chomp;
@f = split(" ", $_);
$pdb_resnum[$f[1]] = $f[4];
}
close F;

print "reading control freqs from $freqfn\n";
$ctrlfreqcount=0;
open FREQF, "$freqfn" || die "cannot open control frequency file $freqfn\n";
while(<FREQF>) {
	chomp;
	#print "read from ctrl file: $_\n";
	($ctrl_mutcode[$ctrlfreqcount], $ctrl_freq[$ctrlfreqcount]) = split(" ", $_);
	#print "parsed into $ctrl_mutcode[$ctrlfreqcount] -- $ctrl_freq[$ctrlfreqcount]\n";
	$ctrlfreqcount++;
	}
close FREQF;
#print "done!\n";

open BIGSCRIPT, "> $tempdir/gnuplot-it-all.sh";

$latexfile="
\\documentclass{article}

\\usepackage{graphicx}

\\oddsidemargin=-0.5in
\\topmargin=-0.5in
\\textwidth=7.8in 
\\textheight=10in 
\\parindent=0pt

\\begin{document}

";

$latexfile = "";

$timept_min = 1000000;
$timept_max = -1000000;

system("grep \"$sample_to_plot \" $tablename > $tablename-tempfile");
$tablename = "$tablename-tempfile";

open F, "$tablename" or die "cannot open $tablename\n";

while(1) {

if (eof(F)) { last; }

# 4 nucleotides times ntimepoints per position
for($ntindex = 0; $ntindex < $numelements; $ntindex++) {

	$ii=0;
	for($i=0;$i<$numtimepointsfromcmdline;$i++) {
		$line  = <F>;
		#print $line;
		chomp($line);

		@temparray = split( " ", $line);
                $samplename = $temparray[0];
                $temptimepoint = $temparray[2];
                $tempmatch = "$samplename-T$temptimepoint";
                $gooddata = 1;
                foreach $x (@blacklist) { if ( $x eq $tempmatch) { $gooddata=0; }}
                #if ($gooddata == 0) { print "match blacklist: $line\n";}   


                if ($gooddata == 1) {
		($samplename, $mutcode[$ii][$ntindex], $timept[$ii][$ntindex], $count[$ii][$ntindex]) = split( " ", $line);
	#	print "$samplename   -- $mutcode[$i][$ntindex] -- $timept[$i][$ntindex] -- $count[$i][$ntindex]\n";

		if ($timept[$ii][$ntindex] > $timept_max) { $timept_max = $timept[$ii][$ntindex]; }
		if ($timept[$ii][$ntindex] < $timept_min) { $timept_min = $timept[$ii][$ntindex]; } 
		$ii++;
		}
		}
	$numactivetimepoints = $ii;
}


print "found $numactivetimepoints active timepoints for $samplename\n"; 
$numtimepoints = $numactivetimepoints;


# figure out timecourse of WT 
for($ntindex = 0; $ntindex < $numelements; $ntindex++) {
	if ($datatype eq "nt") { ($wtnuc, $pos, $mutnt) = $mutcode[0][$ntindex] =~ m/(\D)(\d+)(\D)/; }
	if ($datatype eq "aa") { ($wtnuc, $pos, $mutnt) = $mutcode[0][$ntindex] =~ m/(\D)(\d+)(\D)/; }
	if ($datatype eq "codons") { ($wtnuc, $pos, $mutnt) = $mutcode[0][$ntindex] =~ m/(\D\D\D)(\d+)(\D\D\D)/; }

	if ($mutnt eq $wtnuc) { $wtntindex = $ntindex; last; }	
	}

# calculate the trajectories relative to WT
for($ntindex = 0; $ntindex < $numelements; $ntindex++) {
	for($i=0;$i<$numtimepoints;$i++) {
		# include safeguards against div by zero!
		
		unless ($count[$i][$wtntindex]>0) { print "zero wt count $wtnuc $pos $mutnt  wtindex is $wtntindex\n"; exit(1);}

                if ($count[$i][$ntindex] > 0) {
		$rel_abundance[$i][$ntindex] = log( $count[$i][$ntindex] / $count[$i][$wtntindex] ) /log(2.);
		} else { $rel_abundance[$i][$ntindex] = -30; } # CLAMP AT -30 for -INFINITY # UNFAIR, but should be VISIBLE!
		
	}
}

# create gnuplot scripts for generating every graph

$plotline_cmd = '';

# figure out WT codon and position
($tmp1,$tmp2,$tmp3) = $mutcode[0][0] =~ m/(\D\D\D)(\d+)(\D\D\D)/;

$wtcodon = $tmp1;
$position = $tmp2;


# iterate over all mutant amino acids, in AA order
for($aaidx = 0;$aaidx<length($aalist);$aaidx++) {


	$lscount =0;
	$plotline_cmd = 'plot ';

	$plotctrl_cmd = 'plot ';


	$aaletter = substr($aalist, $aaidx, 1);

	if ($aaletter eq "*") { $aaletter1 = "STOP"; } else {$aaletter1 = $aaletter;}

	$gnufnbase = "$tempdir/plot-$samplename-$position-mut$aaletter1";

# for the amino acid, find all codons and creat data files for plotting

        for($codonidx=0; $codonidx<64; $codonidx++) {
		if ( &getAA($codonlist[$codonidx]) eq $aaletter ) {

			$lscount++;

			$tmp = &getAA($codonlist[$codonidx]);
			#print "codon match: $codonlist[$codonidx] -- $tmp $aaletter\n";		
			# actual codon plot
			$mycodon = $codonlist[$codonidx];
			$mycodonfn = "$tempdir/plot-$position-$mycodon.dat";

			$lscount =1;

			$plotline_cmd .= "\"$mycodonfn\" ls $lscount with linespoints notitle,";

			open DATAFILE, "> $mycodonfn";
			for($ntindex = 0; $ntindex < $numelements; $ntindex++) {
				for($i=0;$i<$numtimepoints;$i++) {
					($wtnuc, $pos, $mutnt) = $mutcode[0][$ntindex] =~ m/(\D\D\D)(\d+)(\D\D\D)/; 
					if ($mutnt eq $mycodon ) { print DATAFILE "$timept[$i][$ntindex] $rel_abundance[$i][$ntindex]\n"; }
					
					}
				}			
			close DATAFILE;	

			# control, noise floor

			$ctrl_F = 0; $lscount = 3;
			for($ii=0;$ii<$ctrlfreqcount;$ii++) {
				$mm  = $wtcodon.$position.$mycodon;
		#		print "control: $mutcode[0][$ntindex]  $ctrl_mutcode[$ii] $ctrl_freq[$ii]\n";
				if ( $ctrl_mutcode[$ii] eq $mm ) { $ctrl_F = $ctrl_freq[$ii]; 
		#				print "got control!\n"; 
						} 
				}
			if ($ctrl_F > 0) { $ctrl_F = log($ctrl_F)/log(2.); } else {$ctrl_F = -30;} # clamp at -30
			$plotctrl_cmd .= "$ctrl_F ls $lscount with lines notitle,";		
			}

	}


	$plotline_cmd =~ s/,$//;
	$plotline_cmd .="\n";

	$plotctrl_cmd =~ s/,$//;
	$plotctrl_cmd .="\n";

        $xrange_min = $timept_min - 0.1*($timept_max - $timept_min);
        $xrange_max = $timept_max + 0.1*($timept_max - $timept_min);

	($tmp1,$tmp2,$tmp3) = $mutcode[0][0] =~ m/(\D\D\D)(\d+)(\D\D\D)/;
#	$mutcode_rewritten = $tmp1."/".&getAA($tmp1).$pdb_resnum[$tmp2].$aaletter;
	$mutcode_rewritten = &getAA($tmp1).$pdb_resnum[$tmp2].$aaletter;

# get doubling time and form x axis tics


	$manualtics = "";
	for($xx = 0; $xx < 10; $xx +=2) {
		$value = "$xx";
		$labelpos = $xx*$doublingtime;
		if ($labelpos < $xrange_max) 
			{$manualtics .= "set label \"$value\" at first $labelpos, -22 font \"Arial,70\"\n"; }
		}

# make one plot for each mutant aa

	open GNUFILE, "> $gnufnbase.gnu";

	print GNUFILE <<END;

	set term post

	set output "$gnufnbase.ps"

	set label \"$samplename $mutcode_rewritten\" at graph 0.1,0.8 font \"Arial,80\"
#	set label \"$samplename $mutcode[0][$ntindex]\" at graph 0.1,0.8 font \"Arial,80\"

	$manualtics

	set xrange[$xrange_min:$xrange_max]
	set yrange[-20:0]

	set xtics font \"Helvetica,50\"
	set ytics font \"Helvetica,50\"

	unset xtics

	set style line 1 lt 1 lw 5 pt 7 ps 5 linecolor rgb \"#ff0000\"
        set style line 2 lt 1 lw 5  linecolor rgb \"#0000ff\"
        set style line 3 lt 1 lw 5  linecolor rgb \"#000000\"

#	set style line 1 lt 1 lw 5  linecolor rgb \"#ff0000\"
#        set style line 2 lt 1 lw 5  linecolor rgb \"#0000ff\"
#        set style line 3 lt 1 lw 5  linecolor rgb \"#00ff00\"
#        set style line 4 lt 1 lw 5  linecolor rgb \"#808000\"
#        set style line 5 lt 1 lw 5  linecolor rgb \"#008080\"
#        set style line 6 lt 1 lw 5  linecolor rgb \"#800080\"


	$plotline_cmd
	$plotfit_cmd
	$plotctrl_cmd

END

	close GNUFILE;

	print BIGSCRIPT  "gnuplot < $gnufnbase.gnu\n";

	$latexfile .= "\\includegraphics[width=3cm,angle=-90]{$gnufnbase.ps}\n";
	$latexfile .= "\\hspace{0.15cm}\n";
	$graphsperline++;

	if ($graphsperline == 4) { $graphsperline =0; $latexfile .= "\n\n\\vspace{5mm}\n\n";}

} # end for ACD...Y*


$latexfile .="\\newpage"; $graphsperline = 0; 

} # end for <F>, reading 64*timepoints blocks

close F;
close BIGSCRIPT;

#$latexfile .= "\\end{document}\n";

open LATEXFILE, "> $plotdir/$outfnbase.tex";
print LATEXFILE $latexfile;
close LATEXFILE;

# Actual printing happens below

print "Making plots\n";

system("sh $tempdir/gnuplot-it-all.sh");
# now we have all the postscript plots

print "Assembling plots into pages\n";

#chdir($tempdir);
#system("latex $outfnbase.tex");
#system("dvips $outfnbase.dvi");
#system("ps2pdf $outfnbase.ps");
#system("cp $outfnbase.pdf $plotdir");

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


