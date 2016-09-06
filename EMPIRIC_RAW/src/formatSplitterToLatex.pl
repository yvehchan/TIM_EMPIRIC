#!/usr/bin/perl

use Cwd;

$outputdir = "../output";
$tempdir ="../temp";
$plotdir ="../plots";
$latexfilename = "splitterlog"; # no extension! (see at the end)

open F, "$outputdir/splitterlog.txt";


while(<F>) {

if ($_ =~ "total reads") { last; }

}

$data = $_;

while(<F>) {

$data .= $_;

}

close F;

#print $data;

$latexfile =  <<LATEX;


\\begin{verbatim}

$data

\\end{verbatim}


LATEX



open F, "> $plotdir/$latexfilename.tex";
print F $latexfile;
close F;

#chdir($tempdir);
#system("latex  $latexfilename.tex");
#system("dvipdf $latexfilename.dvi");
#system("cp     $latexfilename.pdf $plotdir");
