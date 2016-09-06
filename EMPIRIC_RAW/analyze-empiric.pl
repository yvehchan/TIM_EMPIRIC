#!/usr/bin/perl

#print "DID YOU FORGET ulimit -s unlimited so the phred splitter runs ok?";
#sleep(10);

use Cwd;

# For 2013 data
#$fn_fastq = "Jan2014-B1-B8.fastq";
#$fn_fastq = "fake-from-wrongbarcodes.fa" ;
#$fn_barcodes = "barcodes.txt";
#$fn_refseq = "refseq.txt";
#$fn_gene = "SsIGPS-noTag-AA27-248.fa";
#$fn_pdb = "2C3Z.pdb";
#@libraries = ("b1", "b2", "b3","b4","b5", "b6","b7", "b8");
#@rangemin =  (44, 76, 104, 126, 154, 173, 203, 226); # in PDB numbering
#@rangemax =  (53, 85, 113, 135, 163, 182, 212, 235);
#@ctrllibraries = ("b1-wt", "b2-wt", "b3-wt", "b4-wt", "b5-wt", "b6-wt", "b7-wt", "b8-wt", "b1-lib", "b2-lib", "b3-lib", "b4-lib", "b5-lib", "b6-lib", "b7-lib", "b8-lib");  # control/WT MUST have the name
#$inputdir="input";

## Tt-combined - April + August, August barcodes rewritten
if (1)
{
$fn_fastq = "Tt-combined-rewritten-apr-aug-2015.fastq";
$fn_barcodes = "Tt-barcodes-reduced.txt";
$fn_refseq = "Tt-refseq.txt";
$fn_gene = "TtIGPS-noTag-AA-35-254.fa";
$fn_pdb = "1VC4-monomer.pdb";
@libraries = ("b1", "b2", "b3","b4","b5", "b6","b7", "b8");

@rangemin = (44, 78, 106, 128, 153, 174, 207, 229);
@rangemax = (53, 87, 115, 137, 162, 183, 216, 238);

@ctrllibraries = ("b1-wt", "b2-wt", "b3-wt", "b4-wt", "b5-wt", "b6-wt", "b7-wt", "b8-wt", "b1-lib", "b2-lib", "b3-lib", "b4-lib", "b5-lib", "b6-lib", "b7-lib", "b8-lib");

$numtimepoints = 8;
$eattails = 2;  # skip 2 last points from slope calculations
$skiptimezero = "0";

$finalpdf = "empiric-report-Tt-AprAug2015-b1b8.pdf";
$latexfilename = "report-Tt-AprAug2015-b1b8";
$inputdir="input-Tt-combined-Apr-Aug2015";
$extracaption = "Tt/AprAug2015";
$blacklistfn = "Tt-blacklist.txt";
}


# Tt - April 2015
if (0)
{
$fn_fastq = "Apr2015-Tt.fastq";
$fn_barcodes = "Tt-barcodes.txt";
$fn_refseq = "Tt-refseq.txt";
$fn_gene = "TtIGPS-noTag-AA-35-254.fa";
$fn_pdb = "1VC4-monomer.pdb";
@libraries = ("b1", "b2", "b3","b4","b5", "b6","b7", "b8");

@rangemin = (44, 78, 106, 128, 153, 174, 207, 229);
@rangemax = (53, 87, 115, 137, 162, 183, 216, 238);

@ctrllibraries = ("b1-wt", "b2-wt", "b3-wt", "b4-wt", "b5-wt", "b6-wt", "b7-wt", "b8-wt", "b1-lib", "b2-lib", "b3-lib", "b4-lib", "b5-lib", "b6-lib", "b7-lib", "b8-lib");

$numtimepoints = 8;
$eattails = 3;  # skip 2 last points from slope calculations
$skiptimezero = "0";

$finalpdf = "empiric-report-Tt-Apr2015-b1b8.pdf";
$latexfilename = "report-Tt-Apr2015-b1b8";
$inputdir="input-Tt-Apr2015";
$extracaption = "Tt/Apr2015";
$blacklistfn = "Tt-blacklist.txt";
}


# Tm - March 2015
if (0)
{
$fn_fastq = "Tm-mar2015.fastq";
$fn_barcodes = "barcodes.txt";
$fn_refseq = "Tm-refseq.txt";
$fn_gene = "TmIGPS-noTag-AA-32-252.fa";
$fn_pdb = "1I4N-monomer.pdb";
@libraries = ("b1", "b2", "b3","b4","b5", "b6","b7", "b8");
@rangemin =  (40, 74, 102, 124, 150, 172, 202, 224); # in PDB numbering
@rangemax =  (49, 83, 111, 133, 159, 181, 211, 233);
@ctrllibraries = ("b1-wt", "b2-wt", "b3-wt", "b4-wt", "b5-wt", "b6-wt", "b7-wt", "b8-wt", "b1-lib", "b2-lib", "b3-lib", "b4-lib", "b5-lib", "b6-lib", "b7-lib", "b8-lib");
$numtimepoints = 9;
$eattails = 2;  # skip 2 last points from slope calculations
$skiptimezero =0;
$finalpdf = "empiric-report-Tm-b1b8.pdf";
$latexfilename = "report-Tm-b1b8";
$inputdir="input-Tm-mar2015";
$extracaption = "Tm/Mar2015";
$blacklistfn = "Tm-blacklist.txt";
}


if (0)
{
# 2015 data - Ss
$fn_fastq = "L3.fastq";
$fn_barcodes = "barcodes.txt";
$fn_refseq = "Ss-refseq.txt";
#Ss
$fn_gene = "SsIGPS-noTag-AA27-248.fa";
$fn_pdb = "2C3Z.pdb";
@libraries = ("b1", "b2", "b3","b4","b5", "b6","b7", "b8");
@rangemin =  (44, 76, 104, 126, 154, 173, 203, 226); # in PDB numbering
@rangemax =  (53, 85, 113, 135, 163, 182, 212, 235);
@ctrllibraries = ("b1-wt", "b2-wt", "b3-wt", "b4-wt", "b5-wt", "b6-wt", "b7-wt", "b8-wt");
$numtimepoints = 10;
$eattails = 2;  # skip 2 last points from slope calculations
$skiptimezero =1;
$finalpdf = "empiric-report-Ss-b1b8.pdf";
$latexfilename = "report-Ss-b1b8";
$inputdir="input-Ss-jan2015";
$extracaption = "Ss/Jan2015";
}

if (0)
{
# 2015 data - Tt
$fn_fastq = "L4.fastq";
$fn_barcodes = "barcodes.txt";
$fn_refseq = "Tt-refseq.txt";
#Tt

$fn_gene = "TtIGPS-noTag-AA-35-254.fa";  
$fn_pdb = "1VC4-monomer.pdb"; 
@rangemin = (44, 78, 106, 128, 153, 174, 207, 229);
@rangemax = (53, 87, 115, 137, 162, 183, 216, 238);
@libraries = ("b1", "b2", "b3","b4","b5", "b6","b7", "b8");
@ctrllibraries = ("b1-wt", "b2-wt", "b3-wt", "b4-wt", "b5-wt", "b6-wt", "b7-wt", "b8-wt");
$numtimepoints = 10;
$eattails = 2;  # skip 2 last points from slope calculations
$skiptimezero = 1;
$finalpdf = "empiric-report-Tt-b1b8.pdf";
$latexfilename = "report-Tt-b1b8";
$inputdir="input-Tt-jan2015";
$extracaption = "Tt/Jan2015";
}


$do_process_fasta = 1;
$do_create_pdb_mapping = 1;
$do_plot_slopes = 1;
$do_run_gnuplot = 1; # do not re-generate tons of PS files for each aa/codon
$do_heatmap = 1;

$phred_cutoff = 20;


##############################################


$latexfile =  <<LATEX;

\\documentclass{article}
\\usepackage{graphicx}
\\oddsidemargin=-0.5in
\\topmargin=-0.5in
\\textwidth=7.8in 
\\textheight=10in 
\\parindent=0pt
\\begin{document}

{\\LARGE

EMPIRIC Analysis Report

}
LATEX

##############################################

#print "DID YOU FORGET ulimit -s unlimited so the phred splitter runs ok?";
#sleep(10);

chdir("src");


if (1) {


# processing the fasta file
if ($do_process_fasta > 0)
{
$cmd1 = "./processfasta ../$inputdir/$fn_fastq ../$inputdir/$fn_barcodes ../$inputdir/$fn_refseq $phred_cutoff";
print "$cmd1\n";
system($cmd1);
}

#exit(1);

# prepare PDF with processing report (output goes inot plotdir)
system("./formatSplitterToLatex.pl");

$latexfile .= "

\\input{splitterlog.tex}

\\newpage

";

# create position mapping between PDB numbering and libraries
if ($do_create_pdb_mapping > 0) 
{
for($i=0;$i<=$#libraries;$i++) {
	$mylib =   $libraries[$i];
	$myrange1 = $rangemin[$i];
	$myrange2 = $rangemax[$i];
	$cmd2 = "./make-mapping.pl ../$inputdir/$fn_gene ../$inputdir/$fn_pdb $myrange1 $myrange2 $mylib";
	print "$cmd2\n";
	system($cmd2);
	}
}

#exit(1);


# calculate the codon and nt slopes
for($i=0;$i<=$#libraries;$i++) {
	$mylib =   $libraries[$i];
	$cmd3 = "./make-slopes.pl ../output/outtable-codon-bytime.txt $numtimepoints codons $mylib $eattails ../output/maptable-$mylib.txt ../$inputdir/$blacklistfn $skiptimezero"  ;
	print "$cmd3\n";
	system($cmd3);
	$cmd4 = "./make-slopes.pl ../output/outtable-aa-bytime.txt $numtimepoints aa $mylib $eattails ../output/maptable-$mylib.txt ../$inputdir/$blacklistfn $skiptimezero";
	print "$cmd4\n";
	system($cmd4);
	}

# start plotting:
# controls

#exit(1);


for($i=0;$i<=$#ctrllibraries;$i++) {
	$mylib =   $ctrllibraries[$i];
	($myreallib) = $mylib =~ m/(.*)\-.*/;
	$cmd5 = "./plotcontrol.pl ../output/outtable-aa-bytime.txt $mylib aa ../output/maptable-$myreallib.txt";
	system($cmd5);

	$latexfile .="

	\\includegraphics[angle=270,width=6in]{../temp/control-$mylib-aa-plot.ps}

	";


	$cmd5 = "./plotcontrol.pl ../output/outtable-codon-bytime.txt $mylib codons ../output/maptable-$myreallib.txt";
	system($cmd5);

	$latexfile .="

	\\includegraphics[angle=270,width=6in]{../temp/control-$mylib-codons-plot.ps}

	\\newpage
	";

	}

if ($do_plot_slopes > 0) {

# amino acid and codon slopes
for($i=0;$i<=$#libraries;$i++) {
	$mylib =   $libraries[$i];
	$cmd6 = "./plot-slopes.pl ../output/outtable-aa-bytime.txt $numtimepoints aa $mylib ../output/slopes-$mylib-aa.dat ../output/maptable-$mylib.txt $mylib-lib $mylib-wt ../$inputdir/$blacklistfn";
	print "$cmd6\n"; 
	if ( $do_run_gnuplot) {	system($cmd6); }
	$latexfile .= "

	\\input{abundance-$mylib-aa.tex}

	\\newpage
	";

#	$cmd7 = "./plot-slopes.pl ../output/outtable-codon-bytime.txt $numtimepoints codons $mylib ../output/slopes-$mylib-codon.dat ../output/maptable-$mylib.txt $mylib-lib $mylib-wt";
#	if ($do_run_gnuplot ) { system($cmd7); }
#	$latexfile .= "
#
#	\\input{abundance-$mylib-codons.tex}
#
#	\\newpage
#	";


	$cmd8 = "./codonplots.pl ../output/outtable-codon-bytime.txt $numtimepoints codons $mylib ../output/slopes-$mylib-codon.dat ../output/maptable-$mylib.txt $mylib-wt ../$inputdir/$blacklistfn";
	print "$cmd8\n";
	if ($do_run_gnuplot ) { system($cmd8); }

	$latexfile .= "

	\\input{abundance-$mylib-codons-grouped.tex}

	\\newpage
	";

	}

}

}

if ($do_heatmap > 0) {

# generate tables for heatmap and clip the range of fitness values
for($i=0;$i<=$#libraries;$i++) {
	$mylib =   $libraries[$i];
	$cmd9 = "./fitness-reformat.pl ../output/fitness-$mylib-aa.dat ../output/fit-table-$mylib-aa.dat ../output/beta-$mylib.dat";
	system($cmd9);	
	$tableslist .= "../output/fit-table-$mylib-aa.dat ";
	$myliblist .= $mylib;
	$myliblistcomma .= "$mylib,";
	}
	$myliblistcomma =~ s/,$//;
	
#$cmd10 = "./concat-tables.pl $tableslist > ../output/fit-table-$myliblist-aa.dat";
# v2 accepts as many tables as needed..
$cmd10 = "./concat-tables-v2.pl $myliblistcomma  ../output/fit-table-$myliblist-aa.dat";

print "------\n$cmd10\n------\n";

system($cmd10);

# find MmeI sites
$cmd11 = "./find-mmei-sites.pl $myliblistcomma ../$inputdir/$fn_refseq";
system($cmd11);	

# make box file
$cmd12 = "./make-box-file.pl ../output/fit-table-$myliblist-aa.dat $myliblistcomma > ../output/boxfile-$myliblist.dat"; 
system($cmd12); 

# create the labeled heatmap
$cmd13 = "./fitTabletoList.pl ../output/fit-table-$myliblist-aa.dat";
print "$cmd14\n";
system($cmd13);

# create the heatmap without numbers
$cmd14 = "./heatmap-boxes7-nolabels.py  ../output/fit-table-$myliblist-aa.dat ../output/boxfile-$myliblist.dat ../output/fit-table-$myliblist-aa-List.dat ../plots/heatmap-unlabeled-$myliblist.pdf $extracaption";
system($cmd14);

$heatmappdblist .= " heatmap-unlabeled-$myliblist.pdf  ";

$cmd15 = "./heatmap-boxes7.py  ../output/fit-table-$myliblist-aa.dat ../output/boxfile-$myliblist.dat ../output/fit-table-$myliblist-aa-List.dat ../plots/heatmap-labeled-$myliblist.pdf $extracaption";
system($cmd15);

$heatmappdblist .= " heatmap-labeled-$myliblist.pdf  ";

print "$cmd15\n";

# now do heatmaps for each library

for($i=0;$i<=$#libraries;$i++) {

$mylib = $libraries[$i];

$cmd11 = "./find-mmei-sites.pl $mylib ../$inputdir/$fn_refseq";
system($cmd11);

$cmd12 = "./make-box-file.pl ../output/fit-table-$mylib-aa.dat $mylib > ../output/boxfile-$mylib.dat";
system($cmd12);

$onelibrarypdf = "../plots/heatmap-unlabeled-$mylib.pdf";

$cmd14 = "./heatmap-boxes7-onelibrary.py  ../output/fit-table-$mylib-aa.dat ../output/boxfile-$mylib.dat ../output/fit-table-$myliblist-aa-List.dat $onelibrarypdf $extracaption";
system($cmd14);

$allpdfs .= " $onelibrarypdf  ";

$heatmappdblist .= " heatmap-unlabeled-$mylib.pdf  ";

}

$cmd20 = "gs -sDEVICE=pdfwrite -dNOPAUSE -DBATCH -q -sOutputFile=../plots/heatmap-$myliblist-by-library.pdf $allpdfs";
system($cmd20);

#$cmd21 = "./heatmap-boxes7-subplots.py $myliblistcomma ../plots/heatmap-$myliblist-subplots.pdf   $extracaption";
#print "$cmd21\n";
#system($cmd21);
#$heatmappdblist .= "heatmap-$myliblist-subplots.pdf ";


# do not incorporate heatmaps at LaTeX stage
#exit(1);
#$latexfile .="
#
#\\includegraphics[width=7in]{heatmap-$myliblist.ps}
#
#\\newpage
#";


#$latexfile .="
#
#\\includegraphics[width=7in]{heatmap-labeled-$myliblist.ps}
#
#\\newpage
#";

} # end if do_heatmap


# assemble all PDF files together

$pdflist = "splitterlog.pdf ";

for($i=0;$i<=$#ctrllibraries;$i++) {
	$mylib =   $ctrllibraries[$i];
	$pdflist .= "control-$mylib-aa-plot.pdf " ;
	$pdflist .= "control-$mylib-codons-plot.pdf ";
	}

for($i=0;$i<=$#libraries;$i++) {
	$mylib =   $libraries[$i];	
	$pdflist .= "abundance-$mylib-aa.pdf " ;
#	$pdflist .= "abundance-$mylib-codons.pdf ";
	$pdflist .= "abundance-$mylib-codons-grouped.pdf ";
	}

#$pdflist .= "heatmap-$myliblist.pdf ";
#$pdflist .= "heatmap-labeled-$myliblist.pdf ";


$latexfile .= "

\\end{document}

";

print "READY TO ASSEMBLE LATEX FILE\n";

open LF, ">../plots/$latexfilename.tex";
print LF $latexfile;
close LF;


chdir("../plots");
$cmd16 = "latex $latexfilename";
system($cmd16);
system("dvips $latexfilename; ps2pdf $latexfilename.ps");


#exit(1);

$cmd16 = "gs -DBATCH -DNOPAUSE -q -sDEVICE=pdfwrite -sOutputFile=$finalpdf $latexfilename.pdf $heatmappdblist";
print "$cmd16\n";
system($cmd16);

# http://bioeng-student.blogspot.com/2009/02/adding-pdf-headers-footers-and.html
# pdftk report-Ss-b1b8.pdf background ss.pdf output report-Ss-b1b8-watermarked.pdf
