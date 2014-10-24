# Parser to read and quality-filter pyro results

from Bio import SeqIO
from Bio.Seq import Seq
from FastqParser import QualSeq

class Qual:
    def __init__ (self, name, qual):
        self.name = name
        self.quals = qual

    def __getitem__ (self, item):
        return int(quals[item])

    def __str__ (self):
        return "Qual(name = " + self.name + ", quals = [" + (", ".join([str(x) for x in self.quals])) + "]"

class QualFile:
    def __init__ (self, qualfile):
        self.qualfile = qualfile
        try:
            self.line = qualfile.next()
        except StopIteration:
            raise Exception("Format error in qualfile")
        self.eof=0

    def next(self):
        if self.eof == 1:
            raise StopIteration

        quals = [ ]
        
        if not self.line[0] == ">":
            raise Exception("Format error in qualfile")

        seq_name = self.line[1:].split()[0]
        
        while True:
            try:
                self.line = self.qualfile.next()
                
                if self.line[0] == ">":
                    return Qual(seq_name, quals)
        
                quals += [int(x) for x in self.line.split()]
            except StopIteration:
                self.eof = 1
                return(Qual(seq_name, quals))

class RawPyroRes:
    def __init__(self, fasta_file, qual_file=None):
        fasta = open(fasta_file)
        self.qual_present = True
        if qual_file:
            try:
                qual = open(qual_file)
                self.qual = QualFile(qual)
            except IOError:
                print "no qual"
                self.qual_present = False
        else:
            self.qual_present = False
        self.fasta = SeqIO.parse(fasta,"fasta")

    def __iter__ (self):
        return self

    def next(self):
        r = [self.fasta.next(), None]
        seq_record = r[0]
        if self.qual_present:
            qual = self.qual.next()
            if seq_record.id != qual.name:
                raise Exception("Fasta and Qual ID missmatch")

            if len(seq_record.seq) != len(qual.quals):
                raise Exception("Length of sequence and quality mismatch: " + seq_record.id + " " + qual.name)
            r[1] = qual
        return QualSeq(r)






class PyroRes:

    def __init__(self, fasta_file, qual_file, 
                 mean_min=0, min_length=0, 
                 min_qual=0, raw_filtering=0):
        self.rpr = RawPyroRes(fasta_file, qual_file)
        self.check_qual = self.rpr.check_qual
        self.mean_min = int(mean_min)
        self.min_length = int(min_length)
        self.min_qual = int(min_qual)
        self.stats = dict(count = 0,
                          skipped = 0,
                          too_short = 0,
                          low_mean = 0,
                          low_min_quality = 0)
        self.raw_filtering = int(raw_filtering)

    def __iter__ (self):
        return self

    def next (self):
        while True:
            rpr_rec = self.rpr.next()
            seq_record = rpr_rec[0]
            
            self.stats["count"] += 1
            if self.rpr.check_qual:
                qual = rpr_rec[1]

                seq = [ ]
                if self.raw_filtering == 0:
                    seq = filter_full(seq_record, qual, self.min_length, self.mean_min, self.min_qual)
                else:
                    seq = filter_hqr(seq_record, qual, self.min_length, self.mean_min, self.min_qual)
                if seq[0] == None:
                    self.stats[seq[1]] += 1
            #         self.stats["skipped"] += 1
            #         continue

            #     return seq[0]
            # return seq_record
        
