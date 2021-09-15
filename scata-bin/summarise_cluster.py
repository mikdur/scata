#!/usr/bin/env python

import constants, uniseq
import sys, os, re, cPickle
from Bio import AlignIO
from subprocess import call



print "Workdir:", os.getcwd()
print "uname:", " ".join(os.uname())


settings = cPickle.load(open(sys.argv[1] + "/settings.pick"))

tmp_dir = None # os.getenv("TMPDIR")

if tmp_dir == None:
    tmp_dir = settings["work_dir"]

print "Temporary files go in ", tmp_dir

sys.stdout.flush()






clusters = cPickle.load(open(settings["work_dir"] + "/clusters.pick"))

uniseq_to_seq = uniseq.UniseqDB(settings["work_dir"] + "/uniseq_to_seq.pick", "r")

for cl_id in sys.argv[2:]:
    print "Cluster id", cl_id

    cluster = clusters[int(cl_id)]


    if cluster["num"] != int(cl_id):
        raise Exception("Inconsistent data")



    total_size = 0
    rev_count = 0
    genotypes = 0
    singleton_seqs = 0
    reference_seqs = [ ]


    for uniseq in cluster["set"]:
        #print uniseq, uniseq_to_seq[uniseq]["count"]
        total_size += uniseq_to_seq[uniseq]["count"]
        rev_count += uniseq_to_seq[uniseq]["revs"]
        
        for subseq in uniseq_to_seq[uniseq]["seqs"]:
            if "_____ref" in uniseq_to_seq[uniseq]["seqs"][subseq]["tags"]:
                reference_seqs.append({"name" : uniseq_to_seq[uniseq]["seqs"][subseq]["tags"]["_____ref"],
                                       "uniseq" : uniseq,
                                       "subseq" : subseq })
            if not (len(uniseq_to_seq[uniseq]["seqs"][subseq]["tags"].keys()) == 1 and \
                    "_____ref" in uniseq_to_seq[uniseq]["seqs"][subseq]["tags"]):
                genotypes += 1
            if uniseq_to_seq[uniseq]["seqs"][subseq]["count"] == 1:
                singleton_seqs += 1

    homop_factor = genotypes / float(len(list(cluster["set"])))

    proportion_singletons = (singleton_seqs / float(total_size)) if total_size else -1


    # Flatten reference sequence list
    if len(reference_seqs) > 0:
        for r in reference_seqs:
            r["name"] = [a[1] for a in r["name"]]




    # Derive the most common genotypes
    seqs = [ ]
    for uniseq in cluster["set"]:
        for subseq in uniseq_to_seq[uniseq]["seqs"].keys():
            seqs.append([uniseq, subseq, uniseq_to_seq[uniseq]["seqs"][subseq]["count"]])




    seqs.sort(lambda a, b: cmp(b[2], a[2]))
    
    repseq = uniseq_to_seq[seqs[0][0]]["seq"]
    seqs += [["","",0]]* int(settings["num_repseqs"]) # Padding
    seqs = seqs[0:int(settings["num_repseqs"])]



    print "Aligning references and repseqs"
    fas_out = open(tmp_dir + "/repseqs" + cl_id + ".fas", "w")

    cnt_rep = 0
    for (i,a) in enumerate(seqs):
        if a[0]=='':
            continue
        cnt_rep += 1
        fas_out.write(">repseq" + str(i) + " (" + str(a[2]) + " out of " + str(total_size) + ")\n"
                      + a[1] + "\n")

    for r in reference_seqs:
        #print r
        for name in r["name"]:
            fas_out.write(">" + name + " (reference)\n"
                          + r["subseq"] + "\n")
            cnt_rep += 1
    
    fas_out.close()

    repseq_alist = list()
    if cnt_rep > 0:
        # Align the sequences
    
        os.system(settings["muscle"] + " -in " + tmp_dir + "/repseqs" + 
                  cl_id + ".fas -out " + tmp_dir + "/repseqs" + 
                  cl_id + ".fas.aln > " + tmp_dir + "/rs_muscle_out." + cl_id + ".out ") #2>&1")

        alignment = AlignIO.read(open(tmp_dir + "/repseqs" + 
                                      cl_id + ".fas.aln"), "fasta")

    
        repseq_alist = map(lambda a: dict(seq=str(a.seq).upper(), 
                                          id=a.description),
                           list(alignment))
    else:
        print "Empty alignment"
    
    # Flatten reference sequence list
    ref_list = map(lambda a: a["name"], reference_seqs)

    if len(ref_list) > 0:
        ref_list = reduce(lambda a, b: a+b, ref_list)
    refs = ", ".join(ref_list)


    result = dict( total_size = total_size,
                   rev_count = rev_count,
                   proportion_singletons = proportion_singletons, 
                   repseqs = seqs,
                   reference_seqs = reference_seqs,
                   references = refs,
                   repseq_alist = repseq_alist,
                   homop_factor = homop_factor,
                   genotypes = genotypes)

    
    cPickle.dump(result,open(settings["work_dir"] + "/cluster" + cl_id + ".pick",
                        "w"))

    print "Cluster id", cl_id, "done"



