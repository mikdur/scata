#!/usr/bin/env python


import constants, uniseq
import sys, os, pickle


print os.uname()
print sys.argv



settings = pickle.load(open(sys.argv[1] + "/settings.pick"))

outfile = sys.argv[2]


uniseq_to_seq = uniseq.UniseqDB(settings["work_dir"] + "/uniseq_to_seq.pick", "r")

print "Merging clusters..."

clusters = [ ]


for file in sys.argv[3:]:
    print "reading cluster", file
    subcluster = pickle.load(open(file))
    clusters += subcluster


# Collaps clusters to non-redundant clusters


print "Number of clusters: ", len(clusters)
collapsed_clusters = [ ]


while True:
    clusters.sort(lambda a, b: cmp(len(a), len(b)))
    
    print "next", len(clusters), len(collapsed_clusters)
    if len(clusters) == 0:
        break
    c = clusters.pop()

    #print "Size      ", len(c)

    while True:
        to_del = [ ]
        # Walk through clusters and check for overlaps.
        for i in enumerate(clusters):
            if len(c & i[1]):
                #print "overlap is", c & i[1]
                # Ensure to not join on reference
                join = False
                for s in list(c & i[1]):
                    if uniseq_to_seq[s]["count"] > 0:
                        join = True
                    else:
                        #print s, "is a ref, adding to both"
                        c.add(s)
                        i[1].add(s)
                if join:
                    #print "joining"
                    c = c | i[1]
                #else:
                    #print "only refs, not joining"
                to_del.append(i[0])
                #print "overlap", i[0]
        #print "to del", len(to_del)

        to_del.reverse()
        for d in to_del:
            del clusters[d]

        if len(to_del) == 0:
            collapsed_clusters.append(c)
            print "Final size", len(c)
            break
    print "repeating"

print "done"

print "outfile is", outfile

pickle.dump(collapsed_clusters,open(outfile,"w"))
