import sys, time, sge, re, os
from Bio import SeqIO

from constants import *

#motif = "GCATATCAATAAGCGGAGGA"
motif = "TCCTCCGCTTATTGATATGC"
#motif = "AGGAGGCGAATAACTATACG"
bp_to_keep = 1200

errors = ""
warnings = ""

ref_count=0
mean_len=0
max_len=0

n=0
lengths = []
seqs = []
#try:
refs = SeqIO.parse(open(sys.argv[1]),"fasta")
for seq_record in refs:
    n += 1
    spl = re.split(motif,str(seq_record.seq))

    if len(spl) == 1:
        spl = re.split(motif,str(seq_record.seq.reverse_complement()))

    if len(spl) == 1:
        #warnings += "Sequence '%s' (%d, %d bp) does not contain filter motif. Excluded.\n" \
         #           % (seq_record.id, n, len(seq_record.seq))
        continue


    seq = spl[-1][0:bp_to_keep]
    #print seq_record.id, seq

    l=len(seq)
    lengths.append(l)
    seqs.append([seq,seq_record])
    
    if l > ref_max_len:
        errors += "Sequence %d (%s) is too long (%d > %d).\n" \
                  % (n, seq_record.id, l, ref_max_len)
    elif l > ref_long:
        warnings += "Sequence '%s' (%d, %d bp) is longer than %d bp. Consider cutting down.\n" \
                    % (seq_record.id, n, l, ref_long)
    if len(errors) > 10000:
        errors += "Too many errors, bailing out\n"
        break

    if len(warnings) > 20000:
        errors += "Too many warnings. Please revise your reference set.\n"
        break

    if n > 50000:
        errors += "More than 50k references. Is this really a relevant reference set?"
        break


    ref_count = len(lengths)
    mean_len = sum(lengths) / ref_count
    max_len = max(lengths)

        
#except:
#    errors += "Parse error, is it really a FASTA file?\n"

print errors
print "----"
print warnings


print ref_count, mean_len, max_len
