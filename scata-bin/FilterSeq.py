from Bio import SeqIO
from Bio.Seq import Seq


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
