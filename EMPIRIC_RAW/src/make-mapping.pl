#!/usr/bin/perl

use Bio::Seq;
use Bio::SeqIO;


$tempdir = "../temp";

$genefn = $ARGV[0];  # fasta
$pdbfn = $ARGV[1];   # pdf file name

$range_min = $ARGV[2];  # first residue number in pdb file
$range_max = $ARGV[3];  # last residue num

$sample_name = $ARGV[4];

%aa3toaa = (
"ALA", "A",
"ARG", "R",
"ASN", "N",
"ASP", "D",
"CYS", "C",
"GLU", "E",
"GLN", "Q",
"GLY", "G",
"HIS", "H",
"ILE", "I",
"LEU", "L",
"LYS", "K",
"MET", "M",
"PHE", "F",
"PRO", "P",
"SER", "S",
"THR", "T",
"TRP", "W",
"TYR", "Y",
"VAL", "V"
);


my $seqio_object = Bio::SeqIO->new(-file => $ARGV[0]);
my $gene_obj  = $seqio_object->next_seq;

#$gene_obj  = Bio::Seq->new( -seq => $gene, -alphabet => "dna", -display_id => "RECODE_GENE");
$prot_obj = $gene_obj->translate;
$prot_seq_from_fasta = $prot_obj->seq;

print "$prot_seq_from_fasta\n";

#exit(1);

##### reading PDB file (manually, not Bioperl)

$targetchain = "A";

open F, "$pdbfn" or die "cannot open $pdbfn\n";
$i=0; $pdb_seq='';
while(<F>) {
chomp;

#print substr($_,0,4)."\n";

unless ($_ =~ /^ATOM/) { next; }

$atomtype = substr($_, 13,3);
$restype = substr($_,17,3);
$rescode = substr($_,22,5);
$chain = substr($_, 21,1);

$xx = substr($_,30,8);$yy = substr($_,38,8);$zz = substr($_,46,8);

if ($atomtype =~ "CA " && $chain =~ $targetchain ) {
#print "|$atomtype|$restype|$chain|$rescode|$xx|$yy|$zz|\n";

$ResLett[$i]=$aa3toaa{$restype};
$ResType[$i]=$restype;
$ResCode[$i]=$rescode;

$pdb_seq .=$aa3toaa{$restype};

#print "$i $ResCode[$i] $ResType[$i] $ResLett[$i]\n";

$i++;
}

}
$numres_pdb = $i;
close F;

#################################

print "\n$pdb_seq\n";

#################################

# align translated gene and the PDB sequence

$pdb_prot_obj =  Bio::Seq->new( -seq=>$pdb_seq,  -display_id=>"RECODE_PDB");
$gene_prot_obj = Bio::Seq->new( -seq=>$prot_seq_from_fasta,  -display_id=>"RECODE_GENE");


$seqio_protein = Bio::SeqIO->new( -file=>"> $tempdir/tmp-proteins.fa", -format => 'fasta' );
$seqio_protein->write_seq( $pdb_prot_obj );
$seqio_protein->write_seq( $gene_prot_obj);

system("muscle -in $tempdir/tmp-proteins.fa -out $tempdir/tmp-proteins-aligned.fa -maxiters 2");

# read back aligned sequences

$seqio_msa = Bio::SeqIO->new( -file =>"$tempdir/tmp-proteins-aligned.fa", -format =>'fasta');
while( $seq_from_msa = $seqio_msa->next_seq ) {
$disp_id = $seq_from_msa->display_id;
if ($disp_id =~ "RECODE_PDB")  {  $prot_pdb_aligned = $seq_from_msa->seq; } 
if ($disp_id =~ "RECODE_GENE") {  $prot_gene_aligned = $seq_from_msa->seq; }
print "$disp_id\n";
}

print "$prot_pdb_aligned\n$prot_gene_aligned\n";

# pdb_range is known: range_min, range_max

$j = 0; # index over non-gapped aa in pdb sequence
$within_range =0; 
$pdb_seq_cutout = "";
$m = 0; # 0-based index over aa cutout
for($i=0;$i<=length($prot_pdb_aligned); $i++) {

$curr_letter = substr($prot_pdb_aligned,$i,1);

#print "$curr_letter ";

unless ($curr_letter =~ "-") { 
	$current_pdb_resnum = $ResCode[$j] ; 

#	print "$j $ResCode[$j]\n";

	if ($current_pdb_resnum == $range_min) { print "range min $range_min found at i/j $i $j\n"; $range_min_i = $i; $within_range =1; }
	if ($within_range) { $pdb_seq_cutout .= $curr_letter; 
			$cutout_pdb_resnum[$m] = $current_pdb_resnum;
			$cutout_pdb_letter[$m] = $curr_letter;
			$m++;
			}
	if ($current_pdb_resnum == $range_max) { print "range min $range_max found at i/j $i $j\n"; $range_max_i = $i; $within_range =0; }
	}

unless ($curr_letter =~ "-") { $j++; }

}

$j = 0; # index over non-gapped aa in translated gene sequence
$k = 0; # nucleotide index in the gene
$within_range =0;
$gene_cutout = "";
$gene_cutout_aa = "";
$mm = 0;
for($i=0;$i<=length($prot_pdb_aligned); $i++) {

$curr_letter = substr($prot_gene_aligned,$i,1);

unless ($curr_letter =~ "-") { 
     if ($i == $range_min_i) { $within_range = 1; }
     if ($within_range ) { $gene_cutout .= substr( $gene_obj->seq , $k, 3); 
			   $gene_cutout_aa[$mm] .= $curr_letter;
			   $gene_cutout_codon[$mm] = substr( $gene_obj->seq , $k, 3); 
			   $mm++;
			} 
     if ($i == $range_max_i) { $within_range = 0; }
}

unless ($curr_letter =~ "-") { $j++; $k+=3; }

}


# testing
$cutout_obj  = Bio::Seq->new( -seq => $gene_cutout, -alphabet => "dna", -display_id => "CUTOUT");
$prot_cutout_obj = $cutout_obj->translate;
$prot_cutout_testseq = $prot_cutout_obj->seq;

print "gene cut: $gene_cutout\n";
print "  xlated: $prot_cutout_testseq\n";
print "from pdb: $pdb_seq_cutout\n";

$status_ok = 0;
if ($prot_cutout_testseq eq $pdb_seq_cutout) { print "Sequence check OK!\n"; $status_ok=1; }

if ($status_ok) {

	open MAPTABLE, "> ../output/maptable-$sample_name.txt";
	for($i=0;$i<$m;$i++) {
		$cutout_pdb_resnum[$i] =~ s/\s//g; # remove spaces from alt location and other pdb quirks
		print "$sample_name $i $gene_cutout_aa[$i] $gene_cutout_codon[$i] $cutout_pdb_resnum[$i] $cutout_pdb_letter[$i]\n";
		print MAPTABLE "$sample_name $i $gene_cutout_aa[$i] $gene_cutout_codon[$i] $cutout_pdb_resnum[$i] $cutout_pdb_letter[$i]\n";
		}
	close MAPTABLE; 	
	} 
	else { print "mismatches found - cannot write translation table!\n";
	}

