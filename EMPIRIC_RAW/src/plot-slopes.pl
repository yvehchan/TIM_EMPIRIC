#!/usr/bin/perl

use Cwd;
use Switch;

$tablename     = $ARGV[0]; # outtable-XX-by-time.txt
$numtimepointsfromcmdline = $ARGV[1];
$datatype      = $ARGV[2]; # nt | aa | codons
$sample_to_plot = $ARGV[3];
$slopefn        = $ARGV[4]; # slope file for plotting fit data
$mapfn          = $ARGV[5];  # position mapping to PDB file
$ctrlsample1    = $ARGV[6];  # control sample (b3-wt etc) for noise levels
$ctrlsample2    = $ARGV[7];  # control sample (b3-lib etc) for library levels
$blacklistfn    = $ARGV[8]; # blacklist - remove some sample/timepoint 


$aalist = "ILVFYWCMAGNSTHQDEKRP*";

# ver 6: time is in doubling times

$outputdir = "../output";
$plotdir = "../plots"; # for final plots; individual panels are in temp

$plot_fitted_lines = 1;

$outfnbase = "abundance-$sample_to_plot-$datatype";

$freqfn = "$outputdir/freq-$ctrlsample1-$datatype.dat";
$freqfn2 = "$outputdir/freq-$ctrlsample2-$datatype.dat";


$dbltimefn = "$outputdir/dbltime-$sample_to_plot.dat";

$tempdir = "../temp";

if ( $datatype eq "nt" )     { $numelements = 4;}
if ( $datatype eq "aa" )     { $numelements = 21;}
if ( $datatype eq "codons" ) { $numelements = 64;}

open F, "$dbltimefn";
$doublingtime = <F>;
close F;
print "doubling time is $doublingtime\n";


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

print "reading control freqs from $freqfn2\n";
$ctrlfreqcount2=0;
open FREQF, "$freqfn2" || die "cannot open control frequency file $freqfn2\n";
while(<FREQF>) {
	chomp;
	#print "read from ctrl file: $_\n";
	($ctrl_mutcode2[$ctrlfreqcount2], $ctrl_freq2[$ctrlfreqcount2]) = split(" ", $_);
	#print "parsed into $ctrl_mutcode[$ctrlfreqcount] -- $ctrl_freq[$ctrlfreqcount]\n";
	$ctrlfreqcount2++;
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

$latexfile = ""; # preamble is now created by analyze-empiric.pl

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
#		if ($gooddata == 0) { print "match blacklist: $line\n";}
		
		
		if ($gooddata == 1) {
		
		
		($samplename, $mutcode[$ii][$ntindex], $timept[$ii][$ntindex], $count[$ii][$ntindex]) = split( " ", $line);
		#print "$samplename   -- $mutcode[$i][$ntindex] -- $timept[$i][$ntindex] -- $count[$i][$ntindex]\n";
		
		if ($timept[$ii][$ntindex] > $timept_max) { $timept_max = $timept[$ii][$ntindex]; }
		if ($timept[$ii][$ntindex] < $timept_min) { $timept_min = $timept[$ii][$ntindex]; } 
		$ii++;
                }

		}

      $numactivetimepoints = $ii;
}

print "found $numactivetimepoints active timepoints\n";
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
		} else { $rel_abundance[$i][$ntindex] = 0; } # UNFAIR, but should be VISIBLE!
		
	}
}

# create gnuplot scripts for generating every graph

for($ntindex = 0;$ntindex< $numelements; $ntindex++) {

	$fitfn =     "$tempdir/fit-$samplename-$mutcode[0][$ntindex].dat";
	$gnufnbase = "$tempdir/plot-$samplename-$mutcode[0][$ntindex]";

	if ($plot_fitted_lines == 1 ) {
		open SLOPEFILE, "$slopefn";
		while(<SLOPEFILE>) {
			chomp;
			($slope_samp, $slope_mutcode, $slope_slope, $slope_interc, $slope_n, $slope_r2) 
			  	= split(" ", $_);
			if ($slope_samp eq $samplename && $slope_mutcode eq $mutcode[0][$ntindex]) 
				{ #print "$_\n"; 
				last; 
				}
		}
		close SLOPEFILE;

		#print "s,i: $slope_slope $slope_interc\n";

		$xrange_min = $timept_min - 0.1*($timept_max - $timept_min);
		$xrange_max = $timept_max + 0.1*($timept_max - $timept_min);
		                       
		open FITFILE, "> $fitfn";
		for($fit_t = $xrange_min; $fit_t<$xrange_max; $fit_t+=($xrange_max-$xrange_min)/20.) {
			$fit_y = $slope_slope * $fit_t + $slope_interc;
			print FITFILE "$fit_t $fit_y\n";
			}
		close FITFILE;
		$plotfit_cmd = "plot \"$fitfn\" ls 2 with lines notitle\n";
	} else { $plotfit_cmd = ""; } 


        $xrange_min = $timept_min - 0.1*($timept_max - $timept_min);
        $xrange_max = $timept_max + 0.1*($timept_max - $timept_min);

	$yrange_min = -18;
	$yrange_max = -2;
                
	open GNUFILE, "> $gnufnbase.gnu";

	$is_wt_to_wt =0;

	if ($datatype eq "aa") {
		($tmp1,$tmp2,$tmp3) = $mutcode[0][$ntindex] =~ m/(\D)(\d+)(\D)/;
		if ($tmp1 eq $tmp3) { $is_wt_to_wt=1; }
		$mutcode_rewritten = $tmp1.$pdb_resnum[$tmp2].$tmp3;
		#print "labeling: $tmp1,$tmp2,$tmp3 -- $mutcode[0][$ntindex] -- $mutcode_rewritten\n";
		$toplabelsize = 80;
		}
	
	if ($datatype eq "nt") { $mutcode_rewritten = $mutcode[0][$ntindex]; $toplabelsize=80; }

	if ($datatype eq "codons") {
		($tmp1,$tmp2,$tmp3) = $mutcode[0][$ntindex] =~ m/(\D\D\D)(\d+)(\D\D\D)/;
		if ($tmp1 eq $tmp3) { $is_wt_to_wt=1; }
		$mutcode_rewritten = $tmp1.$pdb_resnum[$tmp2].$tmp3;

		$aa1 = &getAA($tmp1);
		$aa2 = &getAA($tmp3);

		$mutcode_rewritten .= " ".$aa1.$pdb_resnum[$tmp2].$aa2;
		$toplabelsize = 45;
		}
	                       

	# figure out frequency from controls
	$ctrl_F = 0;
	for($ii=0;$ii<$ctrlfreqcount;$ii++) {
#		print "control: $mutcode[0][$ntindex]  $ctrl_mutcode[$ii] $ctrl_freq[$ii]\n";
		if ( $ctrl_mutcode[$ii] eq $mutcode[0][$ntindex] ) { $ctrl_F = $ctrl_freq[$ii]; 
#				print "got control!\n"; 
				} 
		}
	if ($ctrl_F > 0) { $ctrl_F = log($ctrl_F)/log(2.); }
	$plotctrl_cmd = "plot $ctrl_F ls 3 with lines notitle\n";	
	$ctrl_F1=$ctrl_F;


	$ctrl_F = 0;
	for($ii=0;$ii<$ctrlfreqcount2;$ii++) {
#		print "control: $mutcode[0][$ntindex]  $ctrl_mutcode[$ii] $ctrl_freq[$ii]\n";
		if ( $ctrl_mutcode2[$ii] eq $mutcode[0][$ntindex] ) { $ctrl_F = $ctrl_freq2[$ii]; 
#				print "got control!\n"; 
				} 
		}
	if ($ctrl_F > 0) { $ctrl_F = log($ctrl_F)/log(2.); }
	$plotctrl_cmd .= "plot $ctrl_F ls 3 with lines notitle\n";	
	$ctrl_F2=$ctrl_F;

	$label_lib_x = $xrange_max*0.2;

	$labeltext_F2 = "noise"; $ctrl_F2new = $ctrl_F2;
	if ($ctrl_F2<$yrange_min) { $ctrl_F2new = $yrange_min+1; $ctrl_F2 = sprintf("%.1f",$ctrl_F2); $labeltext_F2="noise $ctrl_F2 below axis"; } 
	if ($ctrl_F2>$yrange_max) { $ctrl_F2new = $yrange_max-1; 
				    if ($ctrl_F2 != 0) { $ctrl_F2 = sprintf("%.1f",$ctrl_F2); $labeltext_F2="noise $ctrl_F2 above axis";  }
					else { $ctrl_F2new= $yrange_min+1; $labeltext_F2="zero noise";}
					} 
	$ctrl_F2 = $ctrl_F2new;

	$liblabels = "";
	if ($is_wt_to_wt == 0) { $liblabels = "
	set label \"lib\"   at $label_lib_x, $ctrl_F1 font \"Arial,50\" \n
	set label \"$labeltext_F2\" at $label_lib_x, $ctrl_F2 font \"Arial,50\" \n";
	} 

	$manualtics = "";
	for($xx = 0; $xx < 10; $xx +=2) {
		$value = "$xx";
		$labelpos = $xx*$doublingtime;
		if ($labelpos < $xrange_max) 
			{$manualtics .= "set label \"$value\" at first $labelpos, -20 font \"Arial,70\"\n"; }
		}

	print GNUFILE <<END;

	set term post

	set output "$gnufnbase.ps"

	set label \"$samplename $mutcode_rewritten\" at graph 0.1,0.8 font \"Arial,$toplabelsize\"

	$liblabels

	$manualtics

	set xrange[$xrange_min:$xrange_max]
	set yrange[$yrange_min:$yrange_max]

#	set xtics font \"Helvetica,50\"
	set ytics font \"Helvetica,50\"

	unset xtics

	set style line 1 lt 1 lw 5 pt 7 ps 5 linecolor rgb \"#ff0000\"
        set style line 2 lt 1 lw 5  linecolor rgb \"#0000ff\"
        set style line 3 lt 1 lw 5  linecolor rgb \"#000000\"

	$plotfit_cmd
	$plotctrl_cmd

	plot "-" ls 1 with linespoints notitle

END

	for($i=0;$i<$numtimepoints;$i++) {
	if (abs($rel_abundance[$i][$ntindex]) > 1e-10)
	{print GNUFILE "$timept[$i][$ntindex] $rel_abundance[$i][$ntindex]\n"; }
	}

	close GNUFILE;

	print BIGSCRIPT  "gnuplot < $gnufnbase.gnu\n";

	push @graphstoplot, "$samplename-$mutcode[0][$ntindex]";



#	$latexfile .= "\\includegraphics[width=3cm,angle=-90]{$gnufnbase.ps}\n";
#	$latexfile .= "\\hspace{0.15cm}\n";
#	$graphsperline++;
#
#	$graphsperpos++;
#
#	if ($graphsperline == 4) { $graphsperline =0; $latexfile .= "\n\n\\vspace{5mm}\n\n"; }
#
#	if ($datatype eq "codons" && $graphsperpos == 64) { $graphsperpos = 0;  $latexfile .="\\newpage"; }

	}


#if ($datatype eq "aa")  { $latexfile .="\\newpage"; $graphsperline = 0;  } 


}

close F;
close BIGSCRIPT;

# for codon plots, reorder output by mutant aa
@graphstoplot_resorted = ();
if ($datatype eq "codons") {
for($i=0;$i<10;$i++) {
	for($j=0;$j<=length($aalist);$j++) {
	$current_mut_aa = substr($aalist,$j,1);
	foreach $graphcode (@graphstoplot) {
			($code) = $graphcode =~ m/$samplename\-(\D\D\D\d+\D\D\D)/;
			($tmp_wt_codon, $tmp_pos, $tmp_mut_codon) = $code =~ m/(\D\D\D)(\d+)(\D\D\D)/;
			$m_aa = &getAA($tmp_mut_codon);
			#print "code is $code -- $tmp_wt_codon $tmp_pos $tmp_mut_codon $m_aa\n";
			if (  $m_aa eq $current_mut_aa   && $tmp_pos == $i ) { 
				#print "         push $graphcode $code -- $tmp_wt_codon $tmp_pos $tmp_mut_codon $m_aa\n";
				push @graphstoplot_resorted, $graphcode; }
			}
		}
	}
}

# create latex file

$prev_mut_aa = substr($aalist,0,1);
for($i=0;$i<=$#graphstoplot;$i++) {

	if ($datatype eq "aa") { $gnufnbase = "$tempdir/plot-$graphstoplot[$i]"; }
	if ($datatype eq "codons") { $gnufnbase = "$tempdir/plot-$graphstoplot_resorted[$i]"; }

	# for each mutant aa, start a new row of codon plots
	if ($datatype eq "codons") {
			($code) = $graphstoplot_resorted[$i] =~ m/$samplename\-(\D\D\D\d+\D\D\D)/;
			($tmp_wt_codon, $tmp_pos, $tmp_mut_codon) = $code =~ m/(\D\D\D)(\d+)(\D\D\D)/;
			$curr_mut_aa = &getAA($tmp_mut_codon);
			if ($graphsperpos>0  && !($curr_mut_aa eq $prev_mut_aa) ) { 
						$graphsperline=0; 
						$latexfile .= "\n\n\\vspace{5mm}\n\n"; 
						$prev_mut_aa = $curr_mut_aa;
						}
			}


	$latexfile .= "\\includegraphics[width=3cm,angle=-90]{$gnufnbase.ps}\n";
	$latexfile .= "\\hspace{0.15cm}\n";
	$graphsperline++;

	$graphsperpos++;
	if ($graphsperline == 4) { $graphsperline =0; $latexfile .= "\n\n\\vspace{5mm}\n\n"; }


	if ($datatype eq "codons" && $graphsperpos == 64) { $graphsperpos = 0; $graphsperline=0; $latexfile .="\\newpage"; $prev_mut_aa = substr($aalist,0,1); }
	if ($datatype eq "aa" && $graphsperpos == 21)     { $graphsperpos = 0; $graphsperline=0; $latexfile .="\\newpage"; }

}

# latex is now closed from analyze-empiric
#$latexfile .= "\\end{document}\n";

open LATEXFILE, "> $plotdir/$outfnbase.tex";
print LATEXFILE $latexfile;
close LATEXFILE;

# Actual printing happens below

print "Making plots\n";

system("sh $tempdir/gnuplot-it-all.sh");
# now we have all the postscript plots

print "Assembling plots into pages\n";

chdir($tempdir);

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


