#!/usr/bin/perl

$inputfile = $ARGV[0];
$samplename = $ARGV[1];
$datatype = $ARGV[2];  # nt, aa, codon
$mapfn = $ARGV[3]; # mapping to pdb numbering


$outputdir = "../output";
$tempdir = "../temp";
$plotdir = "../plots";

if ($datatype =~ "nt")    { $numelements = 4;  $datalength = 30; }
if ($datatype =~ "aa")    { $numelements = 21; $datalength = 10; }
if ($datatype =~ "codon") { $numelements = 64; $datalength = 10; }


$freqfn = "freq-$samplename-$datatype.dat" ; # frequency file name for megaplotter

$plottype = "freq";
$plotmode = "nowt"; # || haswt


if ($mapfn) {
	#b3 0 S TCT 104 S
	open F, "$mapfn" || die "cannot open $mapfn\n";
	while(<F>) {
		chomp;
		@f = split(" ", $_);
		$pdb_resnum[$f[1]] = $f[4];
		}
	close F;
	}


$i=0;
open F, $inputfile;
while(<F>) {
chomp;
@f = split(" ", $_);

if ($f[0] =~ $samplename) 
{
$mutcode[$i] = $f[1];
$timepoint[$i] = $f[2];
$readcount[$i] = $f[3];
$i++;
}

}
close F;

# calculate frequencies if needed

if ($plottype eq "freq") {

for($k1 = 0; $k1< $datalength; $k1++) {

$wtcount =0;

for($k2 = 0; $k2 < $numelements; $k2++) {
      $idx = $k1*$numelements + $k2;

      if ($datatype =~ "nt" || $datatype =~ "aa" )  {
          ($wt, $pos, $mut) = $mutcode[$idx] =~ m/(\D)(\d+)(\D)/;
          }
      if ($datatype =~ "codon"  )  {
          ($wt, $pos, $mut) = $mutcode[$idx] =~ m/(\D\D\D)(\d+)(\D\D\D)/;
          }
                
      if ($wt eq $mut ) { $wtcount = $readcount[$idx]; }
      }

if ($wtcount>0) {
      for($k2 = 0; $k2 < $numelements; $k2++) {
      $idx = $k1*$numelements + $k2;
      $readcount[$idx] /= $wtcount;
      } 
      }

}

# print out the frequencies in the control for later use in megaplotter
open FREQFILE, "> $outputdir/$freqfn";
for($j=0;$j<$i;$j++) {
	print FREQFILE "$mutcode[$j] $readcount[$j]\n";
	}
close FREQFILE;

} # end if freq


# quash wildtype if needed
if ($plotmode eq "nowt") {
for($j=0;$j<$i;$j++) {
      if ($datatype =~ "nt" || $datatype =~ "aa" )  {
          ($wt, $pos, $mut) = $mutcode[$j] =~ m/(\D)(\d+)(\D)/;
          }
          if ($datatype =~ "codon"  )  {
          ($wt, $pos, $mut) = $mutcode[$j] =~ m/(\D\D\D)(\d+)(\D\D\D)/;
          }
          if ($wt eq $mut) { $readcount[$j] = 0; } 
      }
}

$gnufnbase = "control-$samplename-$datatype-plot";

open GNUFILE, "> $tempdir/$gnufnbase.gnu";

$xmin = -2;
$xmax = $numelements * $datalength + 2;


if ($plottype eq "freq") { $sepbar = 0.1; $label1 = "frequencies"; } else { $sepbar = 1000; $label1="raw counts"; }

$graphlabel = "$label1 -- sample $samplename $datatype";

print GNUFILE <<END;

set term post color enhanced
set output "$tempdir/$gnufnbase.ps"

set xrange[$xmin:$xmax]

set yrange[1e-7:1]
set logsca y

unset xtics

#set label \"$graphlabel\" at graph 0.1,0.9 font \"Arial,20\"

set title \"$graphlabel\"

set style line 1 lt 1 lw 3 pt 7 ps 5 linecolor rgb \"#ff0000\"
set style line 2 lt 1 lw 0.5  linecolor rgb \"#505050\"

set xtics ("" 0) scale 0.0001

END

print GNUFILE  "\n\n";
for($k=0;$k<$datalength*$numelements; $k+= $numelements) {
#print GNUFILE "set label \"$mutcode[$k]\" at $k,-500 rotate by 90\n";

$kk = $k+0.3*$numelements;

if ($datatype =~ "aa"   || $datatype =~ "nt")  {
      ($wt, $pos, $mut) = $mutcode[$k] =~ m/(\D)(\d+)(\D)/;
      }
      
if ($datatype =~ "codon")  {
      ($wt, $pos, $mut) = $mutcode[$k] =~ m/(\D\D\D)(\d+)(\D\D\D)/;
      }
            
if ($datatype =~ "aa" || $datatype =~ "codon") { $pospdb = $pdb_resnum[$pos]; } else { $pospdb = $pos;}

$tick = "$wt$pospdb";
print GNUFILE "set xtics add ( \"$tick\"  $kk) rotate by 90\n";

}

# for amino acids, also print substitutions
if ($datatype eq "aa") {
for($k=0;$k<$datalength*$numelements; $k++) {
if ($readcount[$k]>0) {
$kk = $k-0.5;
print GNUFILE "set label \"$mutcode[$k]\" at $kk,$readcount[$k] rotate by 90 font \"Arial,6\"\n";
}
}
}

print GNUFILE <<END2;

plot "-"  ls 1 with impulses notitle, "-" ls 2 with impulses notitle

END2

for($j=0;$j<$i;$j++) {
print GNUFILE "$readcount[$j]\n"; 
}


print GNUFILE  "e\n";
for($k=0;$k<$datalength*$numelements ; $k+=$numelements ) {
$kk = $k-0.5;
print GNUFILE "$kk $sepbar\n";
}


close GNUFILE;

system("gnuplot < $tempdir/$gnufnbase.gnu"); 
#system("ps2pdf $tempdir/$gnufnbase.ps $plotdir/$gnufnbase.pdf");
