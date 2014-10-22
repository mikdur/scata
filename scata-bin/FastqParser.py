
from Bio import SeqIO
from Bio.Seq import Seq

class Qual:
    def __init__ (self, name, qual):
        self.name = name
        self.quals = qual
    

    def __getitem__ (self, item):
        return int(quals[item])

    def __str__ (self):
        return "Qual(name = " + self.name + ", quals = [" + (", ".join([str(x) for x in self.quals])) + "]"

class Single:
    
    def __init__(self, fastq_file):
        self.fastq = SeqIO.parse(fastq_file, "fastq")
        self.qual_present=True
        

    def __iter__ (self):
        return self

    def next(self):
        rec = self.fastq.next()
        q = Qual(rec.id, rec.letter_annotations["phred_quality"])
        rec.letter_annotations=dict()
        return [rec, q]
        
class Paired:
    
    def __init__(self, fastq_1, fastq_2):
        pass
