#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#include<stdarg.h>

//ver v5a = output reads with barcode mismatches

//ver 7 = output reads with NoMatch reject codes

#define READ_LEN 50

char gencode[]={
//"AAA","AAT","AAG","AAC",
'K', 'N', 'K', 'N',
//"ATA","ATT","ATG","ATC",
'I', 'I', 'M', 'I',
//"AGA","AGT","AGG","AGC",
'R', 'S', 'R', 'S',
//"ACA","ACT","ACG","ACC",
'T', 'T', 'T', 'T',
//"TAA","TAT","TAG","TAC",
   '*', 'Y', '*', 'Y',   // * = STOP
//"TTA","TTT","TTG","TTC",
'L', 'F', 'L', 'F',
//"TGA","TGT","TGG","TGC",
'*', 'C', 'W', 'C',  // * = STOP
//"TCA","TCT","TCG","TCC",
'S', 'S', 'S', 'S',
//"GAA","GAT","GAG","GAC",
'E', 'D', 'E', 'D',
//"GTA","GTT","GTG","GTC",
'V', 'V', 'V', 'V',
//"GGA","GGT","GGG","GGC",
'G', 'G', 'G', 'G',
//"GCA","GCT","GCG","GCC",
'A', 'A', 'A', 'A',
//"CAA","CAT","CAG","CAC",
'Q', 'H', 'Q', 'H',
//"CTA","CTT","CTG","CTC",
'L', 'L', 'L', 'L',
//"CGA","CGT","CGG","CGC",
'R', 'R', 'R', 'R',
//"CCA","CCT","CCG","CCC" 
'P', 'P', 'P', 'P'
};

char codon[][4]={   
"AAA","AAT","AAG","AAC",
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
"CCA","CCT","CCG","CCC"
};

char aalist[] = "ACDEFGHIKLMNPQRSTVWY*";

FILE *fp_textout;

void printfout( char *msg, ...)
{
va_list argptr;

va_start(argptr, msg);
vfprintf(stdout, msg, argptr);
va_end(argptr);
fflush(stdout);

va_start(argptr, msg);
vfprintf(fp_textout, msg, argptr);
va_end(argptr);

}

char translate_codon( char *mycodon) {
int i;
int codon_idx = -1;
for(i=0;i<64;i++) if (!strcmp(mycodon, codon[i])) { codon_idx=i; break; }
if (codon_idx<0) { fprintf(stderr, "bad codon %s\n", mycodon); exit(1); }
return gencode[codon_idx];
}


void translate_sequence( char *aaseq, char *ntseq) {
int i,j=0;
char mycodon[4];
for(i=0;i<strlen(ntseq)/3;i++) {

	mycodon[0]=ntseq[3*i];
	mycodon[1]=ntseq[3*i+1];
	mycodon[2]=ntseq[3*i+2];
	mycodon[3]=0;

	aaseq[j++] = translate_codon(mycodon);
	}
aaseq[j++] = 0;
}



int main(int argc, char **argv)
{
char fn_fastq[100];
char fn_barcodes[100];
char fn_refseq[100];
FILE *fp;

char read_seq[READ_LEN];
char read_phred[READ_LEN];
char read_name[100];
char read_plus[100];

char barcode[100][100];
char uniq_barcode[100][20];
int  timepoint[100];
char sample_name[100][20];
int num_barcodes;
int num_uniq_barcodes;

char refseq[100][100];
char uniq_refseq[100][100];
char refseq_name[100][20];
char refseq_aa[100][100];
int num_refseqs;
int num_uniq_refseqs;

FILE *fastq_output[100][100];

int count_reads =0;
int count_rejects=0;
int reject_code;
int count_reject_by_code[100];
int count_by_uniq_refseq[100];

int read_count_table[50][50];

char *reject_type[] = {"None", "Phred", "NonATGC", "Barcode", "NoMatch", "MultiMatch", "MmeI" };
//			 0        1         2         3          4            5       6

//int snp_table[read_refseq_idx][read_barcode_idx][j][nt_idx]++;
int snp_table[32][20][40][4];
// amino acid count table, up to 32 samples, 20 timepoints, 40 aa long, 21 aa types (20+STOP)
int aa_table[32][20][40][21];
// codon count table
int codon_table[32][20][40][64];

int read_uniq_barcode_idx; // -1 if not a valid barcode, or index in uniq_barcode[][] table
int read_uniq_refseq_idx;

int need_reject; // binary flag if a read has to be rejected

int i,j,k,l,m,s;

//int phred_cutoff = 25;
int phred_cutoff ; //= 28;
int phred_base = 33;  // "magic number" for char to phred decoding

// copy input file names from command line parameters
strcpy(fn_fastq,    argv[1]);
strcpy(fn_barcodes, argv[2]);
strcpy(fn_refseq,   argv[3]);
// read phred cutoff from cmd line
sscanf(argv[4],"%d", &phred_cutoff);


FILE *fp_wrongbar = fopen("../output/wrongbar.txt","w");
FILE *fp_nomatch  = fopen("../output/nomatch.fa", "w"); 
int count_nomatch=1;

// open the log file
fp_textout = fopen("../output/splitterlog.txt","w");

if (!fp_textout) { fprintf(stderr, "cannot write log file\n"); exit(1); }

// read in barcode table;
fp = fopen(fn_barcodes, "r");
if (!fp) { fprintf(stderr, "cannot open %s\n", fn_barcodes); exit(1); }
i=0;
while(!feof(fp)) {
	fscanf(fp, "%s %d %s", barcode[i], &timepoint[i], sample_name[i] );
	i++;	
	}
num_barcodes=i-1;
fclose(fp);

// print the barcode table for test
for(i=0; i < num_barcodes; i++) {
	printfout("idx: %d data: %s %d %s\n", i, barcode[i], timepoint[i], sample_name[i]);
	}

// find unique barcodes
strcpy( uniq_barcode[0], barcode[0]);
num_uniq_barcodes = 1;
for(i=0;i<num_barcodes;i++) {
	int found = 0;
	for(j=0;j<num_uniq_barcodes;j++) 
		if (!strcmp(uniq_barcode[j],barcode[i])) found = 1;
	if (!found) {
		strcpy(uniq_barcode[num_uniq_barcodes], barcode[i]);
		num_uniq_barcodes++;
		}
}

// print the unique barcode table for test
for(i=0; i < num_uniq_barcodes; i++) {
	printfout("idx: %d uniqbarcode: %s \n", i, uniq_barcode[i] );
	}

// read in reference sequence, translate to protein
fp = fopen(fn_refseq, "r");
if (!fp) { fprintf(stderr, "cannot open %s\n", fn_fastq); exit(1); }
i=0;
while(!feof(fp)) {
	fscanf(fp, "%s %s", refseq_name[i], refseq[i] );
	translate_sequence(refseq_aa[i], refseq[i]);
	i++;	
	}
num_refseqs=i-1;
fclose(fp);

// print refernce sequence for tests
printfout("refseq found:\n");
for(i=0; i < num_refseqs; i++) {
	printfout("idx: %d data: %s %s\n", i, refseq_name[i], refseq[i]);
	}

// find uniq reference sequences
strcpy( uniq_refseq[0], refseq[0]);
num_uniq_refseqs = 1;
for(i=0;i<num_refseqs;i++) {
	int found = 0;
	for(j=0;j<num_uniq_refseqs;j++) 
		if (!strcmp(uniq_refseq[j],refseq[i])) found = 1;
	if (!found) {
		strcpy(uniq_refseq[num_uniq_refseqs], refseq[i]);
		num_uniq_refseqs++;
		}
}

printfout("unique refseq:\n");
for(i=0; i < num_uniq_refseqs; i++) {
	printfout("idx: %d data: %s\n", i, uniq_refseq[i]);
	}

// develop mapping between ( uniq_refseq_idx, uniq_barcode_idx ) -> ( samp_name, timepoint )


// open output files;
for(i=0;i<num_uniq_refseqs; i++)
for(j=0;j<num_uniq_barcodes; j++) {
char fastqfilename[100];
sprintf(fastqfilename, "../fastqout/out_refs_%s_barcode_%s.fastq",uniq_refseq[i], uniq_barcode[j]);
fastq_output[i][j] = fopen(fastqfilename, "w");
}


// set arrays to zero

for(i=0;i<100;i++) { 
	count_reject_by_code[i]=0;
	count_by_uniq_refseq[i]=0;
	}

for(i=0;i<32;i++)  
	for(j=0;j<20;j++) 
		for(k=0;k<40;k++) {

		for(l=0;l<4;l++) snp_table[i][j][k][l]=0;
		for(l=0;l<strlen(aalist);l++) aa_table[i][j][k][l]=0;
		for(l=0;l<64;l++) codon_table[i][j][k][l]=0;

		}

for(i=0;i<50;i++) for(j=0;j<50; j++) read_count_table[i][j]=0;

//go through the file

fp = fopen(fn_fastq, "r");
if (!fp) { fprintf(stderr, "cannot open %s\n", fn_fastq); exit(1); }

while(!feof(fp)) {

	if (count_reads %1000000 ==0 ) printf("reads: %.0fM\n", count_reads/1e6);

	fgets(read_name,  100, fp); if (feof(fp)) break;
	fgets(read_seq,   100, fp);
	fgets(read_plus,  100, fp);
	fgets(read_phred, 100, fp);

	read_seq[strlen(read_seq)-1]=0;
	read_phred[strlen(read_phred)-1]=0;
//	read_name[strlen(read_name)-1]=0;
//	read_plus[strlen(read_plus)-1]=0;

	//printf("name : %s\n", read_name);
	//printf("plus : %s\n", read_plus);
	//printf("phred: %s\n", read_phred);

	count_reads++;

// filter by phred scores: 
	need_reject = 0;
	for(i=0;i<strlen(read_phred);i++)
		if (read_phred[i] < phred_base + phred_cutoff) { need_reject=1; break; }

//	need_reject =0;

	if (need_reject) {
		count_rejects++;		
		reject_code = 1;
		count_reject_by_code[reject_code]++;
		// save rejected read?
		continue; // goes back to reading the next short read
		}

// check for any Ns from pos 3 onwards; barcodes can have any letters (for merging datasets)
	for(i=3;i<strlen(read_seq);i++)
		if (read_seq[i]!='A' && read_seq[i]!='T' && read_seq[i]!='G' && read_seq[i]!='C') { need_reject=1; break; }
	if (need_reject) {
		count_rejects++;		
		reject_code = 2;
		count_reject_by_code[reject_code]++;
		continue; // goes back to reading the next short read
		}

// check the barcode
	read_uniq_barcode_idx=-1;
	for(i=0; i < num_uniq_barcodes; i++) {
		
		char read_barcode[10];
		read_barcode[0] = read_seq[0];
		read_barcode[1] = read_seq[1];
		read_barcode[2] = read_seq[2];
		read_barcode[3] = 0;

		if (!strcmp(read_barcode, uniq_barcode[i])) {  read_uniq_barcode_idx=i; break; }

		}

	if (read_uniq_barcode_idx == -1) {
		count_rejects++;		
		reject_code = 3;
		count_reject_by_code[reject_code]++;
		
		fprintf(fp_wrongbar, "%s\n", read_seq);
		
		continue; // goes back to reading the next short read		
		}

// compare read_seq with reference and find the SNP

	int snp_count; // number of SNP of read with respect to current refseq	
	char *read_dna = read_seq + 6; // DNA sequence, after 3 letter barcode and spacer
	char *read_dna_phred = read_phred +6; // Phred scores of read_dna

	int snp_pos[100];
	int match_flag[20];
	
	for( i = 0; i < num_uniq_refseqs; i++ ) {
		snp_count=0;
		match_flag[i] = 0;
		
		for (j = 0; j<30; j++ ) { // scan along the DNA sequence from sample
			if (read_dna[j]!=uniq_refseq[i][j]) { 
							 snp_pos[snp_count] = j;
							 snp_count++;  
							}
			if ( snp_count > 3) break;
 			}
		
		int match_this_refseq = -1 ;

		if (snp_count > 3) match_this_refseq = 0;
		if (snp_count == 3) {  // only ...|NNN|...
				if ( snp_pos[0]%3==0  && snp_pos[1]==snp_pos[0]+1 &&  snp_pos[2]==snp_pos[1]+1 ) 
					{ match_this_refseq = 1; } else { match_this_refseq = 0; }
				}
		if (snp_count == 2) { //  only ...|NN.|... ...|N.N|...  ...|.NN|...
				if ( snp_pos[1]-snp_pos[0] <3 && ( ( snp_pos[0]/3 ) == ( snp_pos[1]/3) ) ) 
					{ match_this_refseq = 1; } else { match_this_refseq = 0; }
				}
		if (snp_count == 1 ) match_this_refseq = 1;
		if (snp_count == 0 ) match_this_refseq = 1;

		if (match_this_refseq == -1 ) { fprintf(stderr,"internal error counting mutants!\n"); exit(1); }
		// here, we know if we match current (i-th) refseq and its EMPIRIC mutants

		match_flag[i] = match_this_refseq;

		}

	int match_count = 0;	
	for(i=0;i<num_uniq_refseqs; i++ ) match_count += match_flag[i];

	if ( match_count == 0 ) {
		count_rejects++;		
		reject_code = 4; //nomatch
		count_reject_by_code[reject_code]++;
		
		count_nomatch++;
		fprintf(fp_nomatch,">r%d\n%s\n",count_nomatch,read_seq);
		
		continue; // goes back to reading the next short read		
		}

	if ( match_count > 1 ) {
		count_rejects++;		
		reject_code = 5; //mutlimatch
		count_reject_by_code[reject_code]++;
		continue; // goes back to reading the next short read		
		}

	read_uniq_refseq_idx = -1;	
	for(i=0;i<num_uniq_refseqs; i++) if ( match_flag[i]==1 ) read_uniq_refseq_idx =i;

	if (read_uniq_refseq_idx == -1 ) { fprintf(stderr, "internal error assigning to refseq\n"); exit(1);}
	
// check for MmeI restriction site and reject TCCRAC + reverse complement
	need_reject = 0;
	char *MmeIsites[] = { "TCCAAC", "TCCGAC",   "GTTGGA", "GTCGGA"};  //  "AGGTTG", "AGGCTG" };
	int mmei_position;
	for( i = 0; i <4 ; i++) {
		char *p = strstr( read_dna, MmeIsites[i]);
		if (p) { need_reject = 1; break; }
		}

	if ( need_reject ) {
		count_rejects++;		
		reject_code = 6;
		count_reject_by_code[reject_code]++;
		continue; // goes back to reading the next short read		
		}

// parsing completed: read_uniq_refseq_idx and read_uniq_barcode_idx are known.

// counting...
	count_by_uniq_refseq[read_uniq_refseq_idx]++;
	read_count_table[read_uniq_refseq_idx][read_uniq_barcode_idx]++;	

// write the read to a separate file by unique bc + unique refseq
	FILE *fpout_fastq = fastq_output[read_uniq_refseq_idx][read_uniq_barcode_idx] ;
	
	fprintf(fpout_fastq, "%s", read_name );
	fprintf(fpout_fastq, "%s\n", read_dna );  // exclude barcode and spacer
	fprintf(fpout_fastq, "%s", read_plus );
	fprintf(fpout_fastq, "%s\n", read_dna_phred ); // matching read_dna


// counting all nucleotides along the DNA	
	for(j=0; j<30; j++) {
		
		int nt_idx;
		if (read_dna[j]=='A') nt_idx = 0;
		if (read_dna[j]=='T') nt_idx = 1;
		if (read_dna[j]=='G') nt_idx = 2;
		if (read_dna[j]=='C') nt_idx = 3;

		// snp_table indices:
		// index of refseq in refseq file
		// index of barcode among unique barcodes, ordered as in barcodes file
		// position (1st nucleotide is position zero) - C string convention
		// nucleotide index (ATGC = 0123)
 		snp_table[read_uniq_refseq_idx][read_uniq_barcode_idx][j][nt_idx]++;

		}


// count all amino acids along the DNA of each read; assume read DNA is in frame

	for(j=0;j<10;j++) {

		// form and translate the codon
		char aa, mycodon[4];
		mycodon[0] = read_dna[3*j];
		mycodon[1] = read_dna[3*j+1];
		mycodon[2] = read_dna[3*j+2];
		mycodon[3] = 0;

		for(k=0;k<64;k++) if (!strcmp(mycodon, codon[k])) break ;
		codon_table[read_uniq_refseq_idx][read_uniq_barcode_idx][j][k]++;

		aa = translate_codon(mycodon);
		for(k=0;k<strlen(aalist);k++) if (aalist[k]==aa) break;

		// k is the amino acid index in ACDEF.. order 		
		aa_table[read_uniq_refseq_idx][read_uniq_barcode_idx][j][k]++;

	}	

};

fclose(fp);

// print statistics:

printfout("total reads: %d, %.1fM\n", count_reads, (double)count_reads/1e6);
printfout("rejected reads: %d, %.1f%%\n", count_rejects, (double)count_rejects/count_reads*100);

printfout("reject reasons:\n");
for(i=1;i<7;i++) printfout("   code %d  %10s, count: %d\n", i, reject_type[i], count_reject_by_code[i]);

printfout("accepted reads by uniq_refseq:\n");
for(i=0;i<num_uniq_refseqs;i++) printfout("   uniq refseq %d, count: %d\n", i, count_by_uniq_refseq[i]);

printfout("accepted reads by sample:\n");

for(i=0;i<num_refseqs;i++) {
	//go over all entries in barcodes file 
	//(barcode-timepoint-refseq combinations)
	for(j=0;j<num_barcodes;j++) {  
		if (!strcmp(refseq_name[i],sample_name[j])) {
			// look up m - index of the sample+timepoint barcode in the uniq barcode array
			m = -1;
			int n;
			for(n=0;n<num_uniq_barcodes;n++) if (!strcmp(uniq_barcode[n],barcode[j])) m = n;
			if (m<0) { fprintf(stderr,"uniq barcode idx - internal error\n"); exit(1); } 

			// look up s - index of the sample+timepoint refseq in the uniq refseq array
			s = -1;
			for(n=0;n<num_uniq_refseqs;n++) if (!strcmp(uniq_refseq[n],refseq[i])) s = n;
			if (s<0) { fprintf(stderr,"unique refseq idx - internal error\n"); exit(1); } 

			printfout( "%s %s time %d, count: %d\n",   
				refseq_name[i], barcode[j], 
				timepoint[j],
				read_count_table[s][m]);
			}
			}

}


// process the snp table

/*
fp = fopen("outtable.txt","w");
for(i=0;i<num_refseqs;i++) {
	//go over all entries in barcodes file 
	for(j=0;j<num_barcodes;j++) {  
		if (!strcmp(refseq_name[i],sample_name[j])) {
			// look up m - index of the current (j-th) barcode in the uniq barcode array
			m = -1;
			for(k=0;k<num_uniq_barcodes;k++) 
				if (!strcmp(uniq_barcode[k],barcode[j])) m = k;
			if (m<0) { fprintf(stderr,"barcode idx - internal error\n"); exit(1); } 
			for(k=0;k<30;k++) {
				fprintf(fp, "%s %s t%d pos%d", refseq_name[i], barcode[j], timepoint[j], k);
				for(l=0;l<4;l++) fprintf(fp," %d", snp_table[i][m][k][l]);
				fprintf(fp,"\n");
				}
			}
		}
	}

fclose(fp);
*/

fp = fopen("../output/outtable-nt-bytime.txt","w");

for(i=0;i<num_refseqs;i++) {
	for(k=0;k<30;k++) {  // 30 nucleotides in a read
		for(l=0;l<4;l++) {
		//go over all entries in barcodes file 
		//(barcode-timepoint-refseq combinations)
			for(j=0;j<num_barcodes;j++) {  
				if (!strcmp(refseq_name[i],sample_name[j])) {

					// look up m - index of the sample+timepoint barcode in the uniq barcode array
					m = -1;
					int n;
					for(n=0;n<num_uniq_barcodes;n++) if (!strcmp(uniq_barcode[n],barcode[j])) m = n;
					if (m<0) { fprintf(stderr,"uniq barcode idx - internal error\n"); exit(1); } 

					// look up s - index of the sample+timepoint refseq in the uniq refseq array
					s = -1;
					for(n=0;n<num_uniq_refseqs;n++) if (!strcmp(uniq_refseq[n],refseq[i])) s = n;
					if (s<0) { fprintf(stderr,"unique refseq idx - internal error\n"); exit(1); } 

					char atgc[] = "ATGC";
					fprintf(fp, "%s %c%d%c %d %d\n",   
						refseq_name[i],
						refseq[i][k], k, atgc[l],
						timepoint[j],
						snp_table[s][m][k][l]);
					}
			}

		}
	}
}

fclose(fp);


fp = fopen("../output/outtable-aa-bytime.txt","w");

for(i=0;i<num_refseqs;i++) {
	for(k=0;k<10;k++) { // 10 amino acids in a read
		for(l=0;l<strlen(aalist);l++) {
		//go over all entries in barcodes file 
		//(barcode-timepoint-refseq combinations)
			for(j=0;j<num_barcodes;j++) {  
				if (!strcmp(refseq_name[i],sample_name[j])) {

					// look up m - index of the sample+timepoint barcode in the uniq barcode array
					m = -1;
					int n;
					for(n=0;n<num_uniq_barcodes;n++) if (!strcmp(uniq_barcode[n],barcode[j])) m = n;
					if (m<0) { fprintf(stderr,"uniq barcode idx - internal error\n"); exit(1); } 

					// look up s - index of the sample+timepoint refseq in the uniq refseq array
					s = -1;
					for(n=0;n<num_uniq_refseqs;n++) if (!strcmp(uniq_refseq[n],refseq[i])) s = n;
					if (s<0) { fprintf(stderr,"unique refseq idx - internal error\n"); exit(1); } 

					//char atgc[] = "ATGC";
					fprintf(fp, "%s %c%d%c %d %d\n",   
						refseq_name[i],
						refseq_aa[i][k], k, aalist[l],
						timepoint[j],
						aa_table[s][m][k][l]);
					}
			}

		}
	}
}


fclose(fp);

fp = fopen("../output/outtable-codon-bytime.txt","w");

for(i=0;i<num_refseqs;i++) {
	for(k=0;k<10;k++) {  // 10 codons in a read
		for(l=0;l<64;l++) { // codon index
		//go over all entries in barcodes file 
		//(barcode-timepoint-refseq combinations)
			for(j=0;j<num_barcodes;j++) {  
				if (!strcmp(refseq_name[i],sample_name[j])) {
					// look up m - index of the current barcode in the uniq barcode array
					m = -1;
					int n;
					for(n=0;n<num_uniq_barcodes;n++) if (!strcmp(uniq_barcode[n],barcode[j])) m = n;
					if (m<0) { fprintf(stderr,"barcode idx - internal error\n"); exit(1); } 

					// look up s - index of the sample+timepoint refseq in the uniq refseq array
					s = -1;
					for(n=0;n<num_uniq_refseqs;n++) if (!strcmp(uniq_refseq[n],refseq[i])) s = n;
					if (s<0) { fprintf(stderr,"unique refseq idx - internal error\n"); exit(1); } 

					char refcodon[4];
					strncpy(refcodon, refseq[i]+3*k, 3);
					refcodon[3] = 0;
					
					fprintf(fp, "%s %s%d%s %d %d\n",   
						refseq_name[i],
						refcodon, k, codon[l],
						timepoint[j],
						codon_table[s][m][k][l]);
					}
			}

		}
	}
}


fclose(fp);

// open output files;
for(i=0;i<num_uniq_refseqs; i++)
for(j=0;j<num_uniq_barcodes; j++) {
fclose(fastq_output[i][j]);
}
 

fclose(fp_wrongbar); 
fclose(fp_nomatch);
fclose(fp_textout);
}

