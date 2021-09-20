#!/usr/bin/env python

import constants, uniseq
import sys, os, re, pickle, subprocess


print(os.uname())
print(sys.argv)


from Bio.Blast import NCBIXML
from Bio.SearchIO._legacy import NCBIStandalone
from subprocess import call

settings = pickle.load(open(sys.argv[1] + "/settings.pick"))

tmp = os.getenv("TMPDIR")

if tmp == None:
    tmp = settings["work_dir"]

print("Temporary files go in ", tmp)

uniseq_to_seq = uniseq.UniseqDB(settings["work_dir"] + "/uniseq_to_seq.pick", "r")

seqs = pickle.load(open(settings["work_dir"] + "/seqs.pick"))

num1 = int(sys.argv[2])
num2 = int(sys.argv[3])
step = int(sys.argv[4])

ids = seqs[num1*step:num1*step+step]
dbids = seqs[num2*step:num2*step+step]

prefix = sys.argv[5]
pick_file = settings["work_dir"] + "/" + prefix

db = tmp + "/" + os.path.basename(prefix) + ".db.fas"
fas = tmp + "/" + os.path.basename(prefix)

print(num1, num2, step)
print(fas)
print(db)
print(pick_file)

dbfile=open(db,"wct")
fasfile=open(fas,"wct")

for id in ids:
    if len(uniseq_to_seq[id]["seq"]) > 10:
        fasfile.write(">" + id + "\n" + uniseq_to_seq[id]["seq"] + "\n" )
for id in dbids:
    if len(uniseq_to_seq[id]["seq"]) > 10:
        dbfile.write(">" + id + "\n" + uniseq_to_seq[id]["seq"] + "\n" )
fasfile.close()
dbfile.close()

formatdb_cmd=subprocess.Popen("which " + settings["formatdb"],shell=True,stdout=subprocess.PIPE).stdout.next()[:-1]
call([formatdb_cmd, "-p", "F", "-i", db])
blastall_cmd=subprocess.Popen("which " + settings["blastall"],shell=True,stdout=subprocess.PIPE).stdout.next()[:-1]

res_h, err_h = NCBIStandalone.blastall(blastall_cmd, "blastn",
                                       db, fas, nprocessors=1,
                                       alignments=500000, 
                                       expectation=settings["blast_expect"])

clusters = { }
cluster_links = list()

blast_records = NCBIXML.parse(res_h)
i = 0
non_hit = ""
analysed_hsps = 0
hit_hsps = 0
for blast_record in blast_records:
    # Take care of non-hit
    if len(non_hit) > 0:
        clusters[non_hit] = set([non_hit])
        

    i += 1
    #print "Doing sequence: ", i
    j=0
    non_hit =  str(blast_record.query)
    for ali in blast_record.alignments:

        for hsp in ali.hsps:
	    analysed_hsps += 1

            query = str(blast_record.query)
            hit = str(ali.hit_def)

            if int("match_from_first" in settings and settings["match_from_first"]) == 1:
                if hsp.query_start != 1 or hsp.sbjct_start !=1: continue
            
            if uniseq_to_seq[query]["count"] == 0 and uniseq_to_seq[hit]["count"] == 0:
                continue

            
            # Fixed length limits
            if float(settings["min_alignment"]) > 1:
                if hsp.align_length < int(settings["min_alignment"]): continue
                if hsp.query_end - hsp.query_start <  int(settings["min_alignment"]): continue
                if hsp.sbjct_end - hsp.sbjct_start <  int(settings["min_alignment"]): continue
            else:
                # Relative hit length limits
                long_len = max(len(uniseq_to_seq[query]["seq"]),
                                len(uniseq_to_seq[hit]["seq"]))
                if len(uniseq_to_seq[query]["seq"]) == long_len  and uniseq_to_seq[query]["count"] > 0:
                    if (float(hsp.query_end) - float(hsp.query_start)) / len(uniseq_to_seq[query]["seq"]) < \
                       float(settings["min_alignment"]): continue
                else:
                    if (float(hsp.sbjct_end) - float( hsp.sbjct_start)) / len(uniseq_to_seq[hit]["seq"]) < \
                       float(settings["min_alignment"]): continue
                    
                
            sbjct_l = list(hsp.sbjct)
            query_l = list(hsp.query)

            # Add divergent sites to distance measure
            distance = sum(map(lambda s, q: float(settings["missmatch_pen"]) if
                                  s != q and (s != "-" and q != "-" and
                                              s != "N" and q != "N") else 0,
                                  sbjct_l, query_l))

            # Skip if distance too long before gap scoring
            if (distance / float(max(hsp.sbjct_end - hsp.sbjct_start,
                                     hsp.query_end - hsp.query_start))
                ) > float(settings["max_dist"]):
                continue

            
            # Add gaps to divergence measure
            # Generator function to generate list of gaps with contexts
            def gaps(a, b, context):
                i = 0
                while i < len(a):
                    if a[i] != "-" and b[i] != "-":
                        i += 1
                        continue
                
                    t = a if a[i] == "-" else b
                    end = i

                    while end < len(a):
                        if t[end] != "-":
                            break
                        end += 1

                    at_end = 1 if i - context < 0 or end + context > len(a) else 0
                    yield [ a[i:end], 
                            a[max(0, i - context):end + context], 
                            b[max(0, i - context):end + context],
                            at_end] \
                            if a[i] == "-" \
                            else [ b[i:end], 
                                   b[max(0, i - context):end + context], 
                                   a[max(0, i - context):end + context],
                                   at_end]
                    i = end
            # Gap scoring function
            def score_gap(gap):
                p = float(settings["gap_extend_pen"])
                op = float(settings["gap_open_pen"])
                hp = int(settings["polyh_num"])
                hw = float(settings["polyh_weight"])
                ew = float(settings["end_gap_weight"])

                score = op + p * len(gap[0])
                if re.match(".*(A{%d}|C{%d}|T{%d}|G{%d}).*" % (hp,hp,hp,hp), 
                            gap[1]) \
                            or re.match(".*(A{%d}|C{%d}|T{%d}|G{%d}).*" % 
                                        (hp,hp,hp,hp),
                                        gap[2]):
                    score = score * hw
                score = score * ew if gap[3] else score
                return score

            # Do the actual scoring
            distance += sum(map(score_gap,
                                gaps(hsp.sbjct,hsp.query,
                                     int(settings["polyh_num"]) - 1)))
             

            distance = distance / float(max(hsp.sbjct_end - hsp.sbjct_start, hsp.query_end - hsp.query_start))
            
            if distance > float(settings["max_dist"]):
                continue
            hit_hsps += 1
            #print "distance", distance, hsp.expect
            j += 1

            # Skip self hit
            if ali.hit_def == blast_record.query: continue

            # Reset non-hit as we have a hit.
            non_hit=""
            
            query = str(blast_record.query)
            hit = str(ali.hit_def)

            query_is_ref_only = False
            hit_is_ref_only = False
            if uniseq_to_seq[query]["count"] == 0:
                query_is_ref_only = True
            if uniseq_to_seq[hit]["count"] == 0:
                hit_is_ref_only = True

            #print i, j, query, hit, query_is_ref_only, hit_is_ref_only

            # Case 1, both query and hit in cluster, this joins
            # two already existing clusters.
            if hit in clusters and query in clusters:


                c1 = clusters[hit]
                c2 = clusters[query]

                if c1 == c2:
                    #print "current - both in same cluster"
                    cluster_links.append(dict( s1 = hit,
                                               s2 = query,
                                               mode = "current",
                                               dist = distance))
                else:
                    if hit_is_ref_only or query_is_ref_only:
                        #print "No join on refs"
                        if not hit_is_ref_only and query_is_ref_only:
                            if hit_is_ref_only:
                                print("ensuring reference hit is in both clusters")
                                if hit not in c1:
                                    c1.add(hit)
                                if hit not in c2:
                                    c2.add(hit)
                            if query_is_ref_only:
                                if query not in c1:
                                    c1.add(query)
                                if query not in c2:
                                    c2.add(query)
                                
                    else:
                        #print "join"
                        # New cluster
                        nc = c1 | c2
                        nc.add(hit)
                        nc.add(query)

                        for id in nc:
                            clusters[id] = nc

                        cluster_links.append(dict( s1 = hit,
                                                   s2 = query,
                                                   mode = "join",
                                                   dist = distance ))

            elif hit in clusters:
                #print "hit"
                if hit_is_ref_only:
                    print("hit is only ref, will not pull in query")
                else:
                    clusters[hit].add(query)
                    clusters[query] = clusters[hit]

                    cluster_links.append(dict( s1 = hit,
                                               s2 = query,
                                               mode = "hit",
                                               dist = distance ))
            elif query in clusters:
                #print "query"
                if query_is_ref_only:
                    print("query is only ref, will not pull in hit")
                else:
                    clusters[query].add(hit)
                    clusters[hit] = clusters[query]
                    cluster_links.append(dict( s1 = hit,
                                               s2 = query,
                                               mode = "query",
                                               dist = distance ))
            else:
                #print "new"

                clusters[hit] = set([hit, query])
                clusters[query] = clusters[hit]
                cluster_links.append(dict( s1 = hit,
                                           s2 = query,
                                           mode = "new",
                                           dist = distance ))

if len(non_hit) > 0:
    clusters[non_hit] = set([non_hit])
    
print("Merging clusters...")


sc = [ ]

for c in list(clusters.values()):
    tt = set()
    for t in sc:
        tt=t
        if c == t:
            break
    if tt == c:
        continue
    sc.append(c)


# Remove zero size clusters (ie cluster with only reference sequences)
to_del = []
for i,c in enumerate(sc):
    sys.stdout.flush()
    if sum([uniseq_to_seq[a]["count"] for a in c]) < 1:
        to_del.append(i)
to_del.reverse()
#print to_del
for d in to_del:
    del sc[d]
    

#map(sum,map(lambda b: map(lambda c: uniseq_to_seq[c]["count"], b), map(lambda a: list(a),c)))

pickle.dump(sc,open(pick_file + ".pick","w"))
if int(settings["graph"]) == 1:
    pickle.dump(cluster_links,open(pick_file + ".links","w"))

print("Analysed HSPs:", analysed_hsps, hit_hsps)
os.remove(db)
os.remove(db + ".nhr")
os.remove(db + ".nsq")
os.remove(db + ".nin")
os.remove(fas)

#print sc
