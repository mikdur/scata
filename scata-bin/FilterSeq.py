from Bio import SeqIO
from Bio.Seq import Seq
from Bio.Alphabet import generic_dna

import difflib


def get_tagset(tagset, tagset_dir):
    if int(tagset) > 0:
        tags = dict()
        tag_file = open( tagset_dir + "/" + str(tagset) + ".txt")
        tag_length = 0
        # Parse tag-file and perform sanity checks
        for line in tag_file:
            line_list = list()
            line = line.rstrip()
            if line == "":
                continue

            line_list = line.split(";")
            #print line_list
            if len(line_list) < 2:
                continue
            tags[line_list[1]] = dict(tag=line_list[0],
                                      matching = set() )
            if tag_length == 0:
                tag_length = len(line_list[1])

            if len(line_list) > 2:
                for i in range(2,len(line_list)):
                    if line_list[i] != "":
                        tags[line_list[1]]['matching'] |= set([line_list[i]])
            
            tags["_____length"] = tag_length
        return tags
    return { }


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

A = 1
C = 2
G = 4
T = 8
trans_table = { 'A' : A,
                     'C' : C,
                     'G' : G,
                     'T' : T,
                     'R' : A | G,
                     'Y' : C | T,
                     'S' : G | C,
                     'W' : A | T,
                     'K' : G | T,
                     'M' : A | C,
                     'B' : C | G | T,
                     'D' : A | G | T,
                     'H' : A | C | T,
                     'V' : A | C | G,
                     'N' : A | C | T | G }


class DeTagSeq:
    # Translation of bases for primer identification.
    def __init__(self, p5, p5s, p3, p3s, t5, t3):
        self.p5s = p5s
        self.p3s = p3s
        self.t5 = t5
        self.t3 = t3
        #self.min_len = min_len

        # Translate primer sequence
        self.p5 = [trans_table[x] for x in p5.upper() if x in trans_table]
        self.p3 = [trans_table[x] for x in p3.upper() if x in trans_table]
        self.p5len = float(len(self.p5))
        self.p3len = float(len(self.p3))


    # Function to detag sequence and return a dict()
    # with sequence and metadata

    def detag_seq(self, seq_record, q=None):
    
        result = dict(rev=False)
        rev = False
        seq = seq_record.seq
        seq_str = str(seq)
        seq_list = [trans_table[x] if not x == 'N' else 0 for x in str(seq).upper() if x in trans_table]
        p5_pos = -1
        accepted_t3 = set()
        if len(self.p5):
            for x in range(0, min(len(seq_str), 2000 + len(self.p5) ) - len(self.p5)):
                m = 0
                for i in range(len(self.p5)):
                    if seq_list[x + i] & self.p5[i]:
                        m += 1
                if m / self.p5len > self.p5s:
                    p5_pos = x
                    break
            
                        
            

            if p5_pos < 0:
                result["rev"] = True
                seq = seq.reverse_complement()
                seq_str = str(seq)
                seq_list = [trans_table[x] if not x == 'N' else 0 for x in str(seq).upper() if x in trans_table]

                for x in range(0, min(len(seq_str), 2000 + len(self.p5) ) - len(self.p5)):
                    m = 0
                    for i in range(len(self.p5)):
                        if seq_list[x + i] & self.p5[i]:
                            m += 1
                    if m / self.p5len > self.p5s:
                        p5_pos = x
                        break
            if p5_pos < 0:
                result["status"] = "no_primer5"
                if q:
                    return (result, None)
                else:
                    return result
        else:
            p5_pos=0
        
        if len(self.t5):
            tag_len = self.t5["_____length"]
            tag_seq = seq_str[(p5_pos - tag_len):p5_pos]
            try:
                result["tag_name"] = self.t5[tag_seq]['tag']
                accepted_t3 = self.t5[tag_seq]['matching']
            except KeyError:
                result["status"] = "no_tag5"
                if q:
                    return (result, None)
                else:
                    return result
        else:
            result["tag_name"] = ""
        

        p3_pos = -1
        p3_matches = [ ]

        if len(self.p3):
            for x in range(len(seq) - len(self.p3) - 1, p5_pos, -1):
                m=0
                for i in range(len(self.p3)):
                    if seq_list[x + i] & self.p3[i]:
                        m += 1
                if m / self.p3len > self.p3s:
                    p3_pos = x
                    break
            if p3_pos < 0:
                result["status"] = "no_primer3"
                if q:
                    return (result, None)
                else:
                    return result
    
            if len(self.t3):
                tag_len = self.t3["_____length"]
                tag_seq = seq_str[(p3_pos + len(self.p3)):(p3_pos + len(self.p3) + tag_len)]
                tag_seq = str(Seq(tag_seq,generic_dna).reverse_complement())
                try:
                    if len(accepted_t3) and self.t3[tag_seq]['tag'] not in accepted_t3:
                        result['status']='chimeric_tag'
                        return (result, None)
                    result["tag_name"] += ("_" + self.t3[tag_seq]['tag'])
                except KeyError:
                    result["status"] = "no_tag3"
                    if q:
                        return (result, None)
                    else:
                        return result
            else:
                result["tag_name"] += ""

            result["seq"] = seq_str[(p5_pos + len(self.p5)):p3_pos]
            seq_record.seq = seq[(p5_pos + len(self.p5)):p3_pos]
            result["seq_record"] = seq_record
        else:
            result["seq"] = seq_str[(p5_pos + len(self.p5)):]
            seq_record.seq = seq[(p5_pos + len(self.p5)):]
            result["seq_record"] = seq_record

        if not q:
            return result
        else:
            if p3_pos:
                q.quals = q.quals[(p5_pos + len(self.p5)):p3_pos]
            else:
                q.quals = q.quals[(p5_pos + len(self.p5)):]
            return (result, q)
# End detag_seq

def old_detag_seq(r, p5, p5s, p3, p3s, t5, t3):
    p5 = p5.upper()
    p3 = p3.upper()
    
    result = dict(rev=False)
    rev = False
    seq = r.seq
    seq_str = str(seq).upper()

    if len(p5):
        p5_matches = [ difflib.SequenceMatcher(None, seq_str[x:x + len(p5)],
                                           p5).ratio() \
                       for x in range(0, min(len(seq_str), 1500) - len(p5)) ]

        if len(p5_matches) == 0:
            result["too_short"] = True
            return result
    

        p5_pos = p5_matches.index(max(p5_matches))

        if p5_matches[p5_pos] < p5s:
            result["rev"] = True
            seq = seq.reverse_complement()
            
            seq_str = str(seq).upper()
            p5_matches = [ difflib.SequenceMatcher(None, seq_str[x:x + len(p5)],
                                               p5).ratio() \
                           for x in range(0, len(seq_str) - len(p5)) ]

            p5_pos = p5_matches.index(max(p5_matches))

        
        if p5_matches[p5_pos] < p5s:
            result["no_5p"] = True
            return result
    else:
        p5_pos=0
        
    if len(t5):
        tag_len = t5["_____length"]
        tag_seq = seq_str[(p5_pos - tag_len):p5_pos]

        try:
            result["tag_name"] = t5[tag_seq]
        except KeyError:
            result["no_such_tag5"] = True
            return result
    else:
        result["tag_name"] = ""
        

    p3_pos = -1
    p3_matches = [ ]

    if len(p3):
        p3_matches = [ difflib.SequenceMatcher(None, seq_str[x:x + len(p3)],
                                               p3).ratio() \
                       for x in range(0,
                                      len(seq_str) - len(p3)) ]

        p3_pos = p3_matches.index(max(p3_matches))

        if p3_matches[p3_pos] < p3s:
            result["no_3p"]=True
            return result

        if len(t3):
            tag_len = t3["_____length"]
            tag_seq = seq_str[(p3_pos + len(p3)):(p3_pos + len(p3) + tag_len)]
            tag_seq = str(Seq(tag_seq,generic_dna).reverse_complement())
            try:
                result["tag_name"] += ("_" + t3[tag_seq])
            except KeyError:
                result["no_such_tag3"] = True
                return result
        else:
            result["tag_name"] += ""

        result["seq"] = seq_str[(p5_pos + len(p5)):p3_pos]
    else:
        result["seq"] = seq_str[(p5_pos + len(p5)):]
        
    return result
# End detag_seq
