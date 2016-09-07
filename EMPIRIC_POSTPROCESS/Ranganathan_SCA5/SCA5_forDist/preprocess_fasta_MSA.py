from Bio import AlignIO

fname = 'ncbi-537-alignment.fasta'
fname_stripped = fname.rstrip('.fasta')


aln = AlignIO.read(fname,'fasta')

N_seq = len(aln)
N_pos = len(aln[0])

gaps_threshold_coeff = 0.8
nonempty_pos = []
for pos in range(N_pos):
    # if num of gaps == N_seq (ALL GAPS), store that index for later removal ...
    # ... DOING THE OPPOSITE, saving nontrivial positions only ...
    # other criteria can be used as well ...
    if (aln[:,pos].count('-') <= gaps_threshold_coeff*N_seq):
        nonempty_pos.append(pos)

# extract nonempty positions from aln and store new aln for manual inspection ...
pos_first = nonempty_pos[0]
extracted_aln = aln[:,pos_first:pos_first+1]
for pos in nonempty_pos[1:]:
    extracted_aln = extracted_aln + aln[:,pos:pos+1] 

# output cleaned fasta ...
AlignIO.write(extracted_aln, fname_stripped+'_cleaned'+'.fasta', "fasta")


# output cleaned free-aln, for MATLAB SCA toolbox ...
with open(fname_stripped+'_cleaned.free','w') as fp:
    for seqrec in extracted_aln:
        fp.write(str(seqrec.seq)+'\n')

