#!/usr/bin/env python

import constants, uniseq
import sys, os, re, pickle, subprocess, time


print(os.uname())
print(sys.argv)


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
print(time.ctime())
for id in ids:
    fasfile.write(">" + id + "\n" + uniseq_to_seq[id]["seq"] + "\n" )
for id in dbids:
    dbfile.write(">" + id + "\n" + uniseq_to_seq[id]["seq"] + "\n" )
fasfile.close()
dbfile.close()

print(time.ctime())

used_ids = dict([(a,0) for a in ids])

usearch_cmd=subprocess.Popen("which " + settings["usearch"],shell=True,stdout=subprocess.PIPE).stdout.next()[:-1]
call([usearch_cmd,
      "--query", fas,
      "--db", db,
      "--evalue", "0.000001",
      "--maxaccepts", "5000",
      "--maxrejects", "5000",
      #"--maxtargets", "50",
      "--id", str(1 - (float(settings["max_dist"]) + 0.1)),
      "--iddef", "0",
      "--userout", (fas + ".out"),
      "--userfields", "query+target+id0+qloz+qhiz+tloz+thiz+ql+tl+opens+exts+pairs+pv+cols",
      "--lopen", str((float(settings["gap_open_pen"]) + 0.0001)),
      "--lext", str((float(settings["gap_extend_pen"]) + 0.0001)),
      "--mismatch", ("-" + str((float(settings["missmatch_pen"]) + 0.0001)))])
      
     



clusters = { }
cluster_links = list()

i = 0
non_hit = ""

col_names = []
analysed_hsps = 0
hit_hsps = 0

res = open(fas + ".out")
for line in res:
    if len(col_names) == 0:
        col_names = line[:-1].split()
        continue

    hit = dict(list(zip(col_names, line[:-1].split())))
    #print col_names,hit    

    query = hit["query"]
    target = hit["target"]
        

    i += 1
    #print "Doing sequence: ", i
    j=0
    analysed_hsps += 1


    if int("match_from_first" in settings and settings["match_from_first"]) == 1:
        if hit["qloz"] != "0" or hit["tloz"] != "0": continue

    # Fixed length limits
    if float(settings["min_alignment"]) > 1:
        if int(hit["cols"]) < int(settings["min_alignment"]): continue

        if int(hit["qhiz"]) - int(hit["qloz"]) <  int(settings["min_alignment"]): continue
        if int(hit["thiz"]) - int(hit["tloz"]) <  int(settings["min_alignment"]): continue
    else: # Relative length limits

        if uniseq_to_seq[query]["count"] == 0 and uniseq_to_seq[target]["count"] == 0:
            continue

        
        long_len = max(len(uniseq_to_seq[query]["seq"]),
                        len(uniseq_to_seq[target]["seq"]))
        if len(uniseq_to_seq[query]["seq"]) == long_len and uniseq_to_seq[query]["count"] > 0:
            if (float(hit["qhiz"]) - float(hit["qloz"])) / len(uniseq_to_seq[query]["seq"]) < \
               float(settings["min_alignment"]): continue
        else:
            if (float(hit["thiz"]) - float(hit["tloz"])) / len(uniseq_to_seq[target]["seq"]) < \
               float(settings["min_alignment"]): continue

    # Add divergent sites to distance measure
    distance = (float(hit["pairs"]) - float(hit["pv"])) * float(settings["missmatch_pen"])

    # Add gaps to distance measure
    distance += float(hit["opens"]) * float(settings["gap_open_pen"])
    distance += float(hit["exts"]) * float(settings["gap_extend_pen"])
    
                    
    distance = distance / float(max(int(hit["thiz"]) - int(hit["tloz"]), int(hit["qhiz"]) - int(hit["qloz"])))

    #print distance, hit["id0"]
            
    if distance > float(settings["max_dist"]):
        continue
    hit_hsps += 1
    #print "distance", distance, hsp.expect
    j += 1


            

    # Skip self hit
    if query == target: continue

    used_ids[query] = 1
    query_is_ref_only = False
    target_is_ref_only = False
    
    if uniseq_to_seq[query]["count"] == 0:
        query_is_ref_only = True
    if uniseq_to_seq[target]["count"] == 0:
        target_is_ref_only = True

    #print i, j, query, target, query_is_ref_only, target_is_ref_only

    # Case 1, both query and target in cluster, this joins
    # two already existing clusters.
    if target in clusters and query in clusters:


        c1 = clusters[target]
        c2 = clusters[query]

        if c1 == c2:
            #print "current - both in same cluster"
            cluster_links.append(dict( s1 = target,
            s2 = query,
            mode = "current",
            dist = distance))
        else:
            if target_is_ref_only or query_is_ref_only:
                #print "No join on refs"
                if not target_is_ref_only and query_is_ref_only:
                    if target_is_ref_only:
                        print("ensuring reference target is in both clusters")
                        if target not in c1:
                            c1.add(target)
                        if target not in c2:
                            c2.add(target)
                    if query_is_ref_only:
                        if query not in c1:
                            c1.add(query)
                        if query not in c2:
                            c2.add(query)
                                
                else:
                    #print "join"
                    # New cluster
                    nc = c1 | c2
                    nc.add(target)
                    nc.add(query)

                    for id in nc:
                        clusters[id] = nc


    elif target in clusters:
        #print "target"
        if target_is_ref_only:
            print("target is only ref, will not pull in query")
        else:
            clusters[target].add(query)
            clusters[query] = clusters[target]
            
            cluster_links.append(dict( s1 = target,
                                       s2 = query,
                                       mode = "target",
                                       dist = distance ))
    elif query in clusters:
        #print "query"
        if query_is_ref_only:
            print("query is only ref, will not pull in target")
        else:
            clusters[query].add(target)
            clusters[target] = clusters[query]
            cluster_links.append(dict( s1 = target,
                                       s2 = query,
                                       mode = "query",
                                       dist = distance ))
    else:
        #print "new"
        
        clusters[target] = set([target, query])
        clusters[query] = clusters[target]
        cluster_links.append(dict( s1 = target,
                                   s2 = query,
                                   mode = "new",
                                   dist = distance ))


    
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
    if sum([uniseq_to_seq[a]["count"] for a in c]) < 1:
        to_del.append(i)
to_del.reverse()
#print to_del
for d in to_del:
    del sc[d]

    
singletons = [ ]
for id,used in used_ids.items():
    if used == 0:
        if uniseq_to_seq[id]["count"] != 0:
            singletons.append(set([id]))


#print singletons

sc += singletons


pickle.dump(sc,open(pick_file + ".pick","w"))
if int(settings["graph"]) == 1:
    pickle.dump(cluster_links,open(pick_file + ".links","w"))

print("Analysed HSPs:", analysed_hsps, hit_hsps)
#os.remove(db)
#os.remove(fas + ".out")
#os.remove(fas)

#print sc
