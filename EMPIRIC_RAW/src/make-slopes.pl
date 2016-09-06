#!/usr/bin/perl

$tablename = $ARGV[0];
$numtimepointsfromcmdline = $ARGV[1];
$datatype = $ARGV[2]; # nt | aa | codons
$sample_to_plot = $ARGV[3];
$eattails = $ARGV[4];
$mapfn = $ARGV[5];  # REQUIRED for amino acid!
$blacklistfn = $ARGV[6]; 
$skiptimezero = $ARGV[7]; # 0/1 to skip first time point from slope calculation


$outputdir = "../output";

# example: ./make-slopes.pl ../outtable-aa-bytime.txt 4 aa b3

# Algo details:
# Skip timepoint t = 0 
# eattails only "eaten" if 3+ points remain.

$tempdir = "tmp";

if ( $datatype eq "nt" )     { $numelements = 4;  $numpositions = 30;}
if ( $datatype eq "aa" )     { $numelements = 21; $numpositions = 10;}
if ( $datatype eq "codons" ) { $numelements = 64; $numpositions = 10; }

$slopefn = "slopes-$sample_to_plot-$datatype.dat";
$slopewarnfn = "slopes-warn-$sample_to_plot-$datatype.dat";
$normslopefn = "fitness-$sample_to_plot-$datatype.dat";
$normcoefffn = "slopeSTOP-$sample_to_plot.dat";  # generated in codons mode
$dbltimefn = "dbltime-$sample_to_plot.dat" ; # has doubling time


$blacklistlen=0;
open F, "$blacklistfn";
while(<F>) {
chomp;
push @blacklist, $_;
}
close F;
foreach $x (@blacklist) { print "BLACKLISTED: $x\n"; }

if (! $mapfn && $datatype eq "aa" ) { die "must provide mapping file for amino acids!\n"; }

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


open SLOPEFILE, "> $outputdir/$slopefn" or die "cannot write slope file\n";
open SLOPEWARN, "> $outputdir/$slopewarnfn" or die "cannot write slope file\n";

system("grep \"$sample_to_plot \" $tablename > $tablename-tempfile");
$tablename = "$tablename-tempfile"; 

open F, "$tablename" or die "cannot open $tablename\n";


while(1) {

if (eof(F)) { last; }

# 4($numelements) nucleotides times ntimepoints per position
for($ntindex = 0; $ntindex < $numelements; $ntindex++) {

        $ii = 0;
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
		$ii++;
		}
               
		#print "$samplename   -- $mutcode[$i][$ntindex] -- $timept[$i][$ntindex] -- $count[$i][$ntindex]\n";
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

		$rel_abundance[$i][$ntindex] = $count[$i][$ntindex] / $count[$i][$wtntindex];
	}
}

#####################################################
#
# $rel_abundance[$i][$ntindex]  $i = 0.. $numtimepoints, $ntindex = 0..3 / 0..20 / 0..63
#
######################################################

# calculate slopes

@slopes=();


for($ntindex = 0;$ntindex<$numelements; $ntindex++) {

	print "$mutcode[0][$ntindex]  ";

	$slope[$ntindex] = 0;
	$intercept[$ntindex] = 0;
	$r2[$ntindex]=0;
	

	# create arrays with x and y for slope determination

	@xvalue=();
	@yvalue=();

	for($i=0;$i<$numtimepoints;$i++) {

		#print "$timept[$i][$ntindex] $rel_abundance[$i][$ntindex]\n";

	if ($rel_abundance[$i][$ntindex] > 0 ) {
	             # && $timept[$i][$ntindex] > 0) 
	        unless ( $timept[$i][$ntindex] ==0 && $skiptimezero eq "1" )      
                {  # SKIP TIME 0!
		push( @xvalue, $timept[$i][$ntindex]);
		push( @yvalue, log($rel_abundance[$i][$ntindex])/log(2.));
		}
		}
	
	}

	$numpushed = $#xvalue+1;

	print "x,y values ready - $numpushed pts  ";

	if ($numtimepoints > $numpushed) { print "warning: zero abundance skipped, $numtimepoints -> $numpushed\n"; }

	if ($numpushed > 1) {	

		($slope[$ntindex], $intercept[$ntindex], $r2[$ntindex]) = &linear_fit( \@xvalue, \@yvalue, \$r2 ); 
			
		} else { print "error: not enough data points $numpushed for this subst\n"; } 

	print "slope = $slope[$ntindex]\n";
	

	$zeros_count = 0;

	for($i=0;$i<$numtimepoints;$i++) {
		if ($rel_abundance[$i][$ntindex] == 0 ) { $zeros_count++; }
		}	
	if ($zeros_count>2) {
			if ($datatype eq "aa" ) {
			($tmpwt, $tmppos, $tmpmut) = $mutcode[0][$ntindex] =~ m/(\D)(\d+)(\D)/;
			$tmppos = $pdb_resnum[$tmppos];
			$tmppos = $tmpwt.$tmppos.$tmpmut;
			print SLOPEWARN "$samplename $tmppos\n";
			}			
			}

	print SLOPEFILE "$samplename $mutcode[0][$ntindex] $slope[$ntindex] $intercept[$ntindex] $numpushed $r2[$ntindex]\n";
}




}

close F;
close SLOPEWARN;
close SLOPEFILE;


if ($datatype eq "codons") {

$avg_slope_of_STOP = 0;
$avg_s2=0;
$stop_count = 0;

@stopslope_array = ();

	open F, "$outputdir/$slopefn";
	for($j=0; $j<$numpositions; $j++){
		for($i=0; $i<$numelements;$i++) {

			$line[$i] = <F>;
			@f = split(" " , $line[$i]);
			($wt, $pos, $mut) = $f[1] =~ m/(\D\D\D)(\d+)(\D\D\D)/;
			if ($mut eq "TAA" ||  $mut eq "TAG" || $mut eq "TGA") { 
				print $line[$i]; 
				#if ($f[2]<0) 
					{ # assume positive slopes of STOP are artifacts..
					push (@stopslope_array, $f[2]);
					$avg_s2 += ($f[2]*$f[2]);
					$stop_count++; 
					$avg_slope_of_STOP += $f[2]; 
					}
				}
				
		}
	}

	close F;

	$avg_slope_of_STOP /= $stop_count;
	$stdev_of_STOP = sqrt( $avg_s2/$stop_count -  $avg_slope_of_STOP*$avg_slope_of_STOP );
	print "Average slope of STOP: $avg_slope_of_STOP, stdev $stdev_of_STOP \n";
	$dbl_time = 1./ $avg_slope_of_STOP;
	print "Doubling time ".sprintf("%.3f\n", (-1)*$dbl_time);

	@stopslope_sorted = sort {$a <=> $b} @stopslope_array;
	$num_elem = $#stopslope_sorted;
	$num_half = $num_elem/2;
	$median = $stopslope_sorted[$num_half];	

	print "MEDIAN slope of STOP: $num_elem $median\n";

	for($iiii=0; $iiii<=$num_elem;$iiii++) { print "SORT $iiii $stopslope_sorted[$iiii]\n"; }


	open DBLFILE, "> $outputdir/$dbltimefn";
	print DBLFILE sprintf("%.3f\n", (-1)*$dbl_time);
	close DBLFILE;
	

	open NORMCOEFF, "> $outputdir/$normcoefffn";
	print NORMCOEFF "$avg_slope_of_STOP\n";
	close NORMCOEFF;

}

### re-normalizing slopes to -1 at missense

if ($datatype eq "aa") {

	open NORMCOEFF, "$outputdir/$normcoefffn" || die "No STOP codon slope file found. Please process the codon table first!\n";
	$missenseslope = <NORMCOEFF>;
	close NORMCOEFF;	

if (abs($missenseslope) < 1e-10) { print "ERROR: MISSENSE SLOPE ZERO ".sprintf("%f", $missenseslope)."\n"; exit(1);}

	open F, "$outputdir/$slopefn";
	open NORMSLOPE, "> $outputdir/$normslopefn";
	for($j=0; $j<$numpositions; $j++){

		for($i=0; $i<$numelements;$i++) {
			$line[$i] = <F>;
			@f = split(" " , $line[$i]);
			$f[2] /= (-$missenseslope);
			($wt, $pos, $mut) = $f[1] =~ m/(\D)(\d+)(\D)/;
			$pos = $pdb_resnum[$pos];
			$f[1] = $wt.$pos.$mut;
			print NORMSLOPE "$f[0] $f[1] $f[2]\n";
			}

	}
	close F;
	close NORMSLOPE;
}



sub linear_fit() {

my @x = @{$_[0]};
my @y = @{$_[1]};
my $i;

my $n = $#x+1;

if ($n>$eattails+2) { $n -= $eattails; }

print "array size is $n\n";

my $xavg=0;
my $yavg=0;
my $xyavg=0;
my $x2avg=0;

for($i=0;$i<$n;$i++)
{
$xavg+=$x[$i];
$yavg+=$y[$i];
$xyavg += ($x[$i]*$y[$i]);
$x2avg += ($x[$i]*$x[$i]);
}

$xavg/=$n;
$yavg/=$n;
$xyavg/=$n;
$x2avg/=$n;

my $ssxx = 0;
my $ssyy = 0;
my $ssxy = 0;

for($i=0;$i<$n;$i++)
{
$ssxx += ($x[$i]-$xavg)*($x[$i]-$xavg);
$ssyy += ($y[$i]-$yavg)*($y[$i]-$yavg);
$ssxy += ($x[$i]-$xavg)*($y[$i]-$yavg);
}

my $r2;
if ($ssyy>0 && $ssxx>0) {
$r2 = ($ssxy*$ssxy/($ssxx*$ssyy));
} else {$r2=0;}

#print "--- $r2\n";

my $slope = ($xyavg -$xavg*$yavg)/( $x2avg- $xavg*$xavg);
my $intercept = $yavg - $slope * $xavg;

return ($slope, $intercept, $r2);

}


