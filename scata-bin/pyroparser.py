# Parser to read and quality-filter pyro results

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
    def __init__(self, fasta_file, qual_file):
        fasta = open(fasta_file)
        self.check_qual = True
        try:
            qual = open(qual_file)
            self.qual = QualFile(qual)
        except IOError:
            print "no qual"
            self.check_qual = False

        self.fasta = SeqIO.parse(fasta,"fasta")

    def __iter__ (self):
        return self

    def next(self):
        r = [self.fasta.next(), None]
	seq_record = r[0]
        if self.check_qual:
            qual = self.qual.next()
            if seq_record.id != qual.name:
                raise Exception("FNA and Qual missmatch")

            if len(seq_record.seq) != len(qual.quals):
                raise Exception("Length of sequence and quality mismatch: " + seq_record.id + " " + qual.name)
            r[1] = qual
        return r

def filter_full(seq_record, qual, min_length, mean_min, min_qual):
    if len(seq_record.seq) < min_length:
        return [None, "too_short"]
        
    if float(sum(qual.quals)) / len(qual.quals) < mean_min:
        return [None, "low_mean_quality"]
            
    seq_list = list(seq_record.seq)
            
    seq_list = map(lambda b, q: b if q >= min_qual else 'N',
                   seq_list, qual.quals)
    
    if seq_list.count('N'):
        return [None, "low_min_quality"]
    
    seq_record.seq = Seq("".join(seq_list), seq_record.seq.alphabet)
    
    return [seq_record, None]


def filter_hqr(seq_record, qual, min_length, mean_min, min_qual):
    if len(seq_record.seq) < min_length:
        return [None, "too_short"]

    # Get regions with quality above minimum
                    
    hqr = [ ]
                    
    t_start=-1
    t_end=-1

    for i in range(len(qual.quals)):
        if t_start == -1: # Outside region
            if qual.quals[i] >= min_qual:
                t_start = i
        else: # Inside region
            if qual.quals[i] < min_qual:
                t_end = i
                hqr.append(dict(start = t_start,
                                end = t_end,
                                mean_quality = 0,
                                length = 0))
                t_start = -1
                t_end = -1
    if t_start != -1 and t_end == -1:
        hqr.append(dict(start = t_start,
                        end = len(qual.quals) - 1,
                        mean_quality = 0,
                        length = 0))
    # Check and modify regions to ensure mean quality is above limit
    #print "Before:", hqr
    for r in hqr:
        if r["end"] - r["start"] < min_length:
            continue
        while r["end"] > r["start"] and \
                (sum(qual.quals[r["start"]:r["end"]]) / float(r["end"] - r["start"])) \
                < mean_min:
            if qual.quals[r["start"]] < qual.quals[r["end"]]:
                r["start"] += 1
            else:
                r["end"] -= 1
        if r["end"] > r["start"]:
            r["mean_qual"] = (sum(qual.quals[r["start"]:r["end"]]) / float(r["end"] - r["start"]))
            r["length"] = r["end"] - r["start"]
    #print "After:", hqr

    hqr = filter(lambda a: a["length"] >= min_length and a["mean_qual"] >= mean_min,
                             hqr)
    
    if not len(hqr):
        return [None, "low_mean_quality"]

    hqr.sort(key=lambda a: a["length"], reverse=True)

                    
    if hqr[0]["end"] - hqr[0]["start"] < min_length:
        return [None, "low_mean_quality"]
                    
    seq_record.seq = seq_record.seq[hqr[0]["start"]:hqr[0]["end"]]
    return [seq_record, None]





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
                    self.stats["skipped"] += 1
                    continue

                return seq[0]
            return seq_record
        
