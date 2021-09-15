#!/usr/bin/env python
import constants
from uniseq import UniseqDB
import sys, os, shutil, math
from subprocess import call
from PyroParser import PyroRes
from sge import SGEJob as GridJob
import cPickle
import re, random
from Bio import SeqIO





def reduce_homopolymer(seq, n):
    if n:
        seq = re.sub("A{%d,}" % (n), "A" * n, seq)
        seq = re.sub("C{%d,}" % (n), "C" * n, seq)
        seq = re.sub("T{%d,}" % (n), "T" * n, seq)
        seq = re.sub("G{%d,}" % (n), "G" * n, seq)    
    return seq

def run_scata(config_file,pr=None):
    # Open settings file and derive all run parameters from it
    settings = { }
    settings["settings_file"] = config_file
    set_in = open(settings["settings_file"])
    for line in set_in:
        tmp = line.split("\n")[0].split(": ")
        if len(tmp) == 2:
            settings[tmp[0]] = tmp[1]
    set_in.close()


    # Create workingdir now, to be able to recover from interruptions.
    os.makedirs(settings["work_dir"])
    os.makedirs(settings["output_dir"])
    os.makedirs(settings["output_dir"] + "/clusters")
    os.makedirs(settings["output_dir"] + "/repseq_alignments")
    os.makedirs(settings["output_dir"] + "/tag_summaries")
    os.makedirs(settings["output_dir"] + "/tag_fastafiles")
    os.makedirs(settings["output_dir"] + "/tag_mappings")

    cPickle.dump(settings,open(settings["work_dir"] + "/settings.pick", "w"))

    cpu = 0.0
    wallclock = 0.0




    print "Reading pyro results and doing initial redundancy screening..."
    seqs = { }
    used_tags = dict()
    seq_stats = { }
    global_seqs = 0

    for seq_a in zip( settings["454names"].split(),
                    zip(settings["454seq"].split(), settings["454stat"].split())):
        # Load datasets
        t_seqs = cPickle.load(open(seq_a[1][0]))
        print "Adding ", seq_a[0]
        added_seqs = 0
        for (k, v) in t_seqs.iteritems():
            k = k.upper()
            global_seqs += len(v)
            added_seqs += len(v)
            seqs[k] = seqs.get(k, []) + v
            for t in v:
                used_tags[t[0]] = used_tags.get(t[0],0) + 1
        print "added %d sequences" % ( added_seqs )
        
        # Load statistics
        seq_stats[seq_a[0]] = cPickle.load(open(seq_a[1][1]))


    print "Number of unique genotypes: ", len(seqs)
    if len(seqs) == 0:
        cluster_summary = open(settings["output_dir"] + "/all_clusters.txt","w")

        cluster_summary.write("General statistics for clustering run " +
                          settings["job_name"] + "\n")
        cluster_summary.write("No sequences to cluster.")

    if settings.get("remove_lowfreq", 0):
        print "Removing low frequency genotypes (< %s)" % settings["remove_lowfreq"]
        new_seqs = dict()
        lf = int(settings["remove_lowfreq"])
        for k, v in seqs.iteritems():
            if len(v) > lf:
                new_seqs[k] = v
        seqs = new_seqs
        new_seqs = dict()
        print "Number of unique genotypes after low frequency pruning: ", len(seqs)

    print "Dumping tag/seqID mappings"
    tag_mapping = { }
    for (k, v) in seqs.iteritems():
        for s in v:
            if len(s) >= 4:
                if s[3][0] not in tag_mapping:
                    tag_mapping[s[3][0]] = { s[3][1] : [s[3][2]] }
                else:
                    tag_mapping[s[3][0]][s[3][1]] = tag_mapping[s[3][0]].get(s[3][1], []) + [ s[3][2] ]

    if len(tag_mapping):
        for (d, v) in tag_mapping.iteritems():
            for (t, s) in v.iteritems():
                f = open(settings["output_dir"] + "/tag_mappings/" + d.replace(" ","_").replace("/","_") + "-" + 
                     t.replace(" ","_").replace("/","_") + ".txt", "w")
                for s_id in s:
                    f.write(s_id + "\n")
                f.close()

    # Downsample tags to given size
    if len(tag_mapping) and int(settings.get("downsample_size", 0)) > 0:
        print "Downsampling tags to at most %s reads" % settings["downsample_size"]

        # Generate subsets
        for d, v in tag_mapping.iteritems():
            for tag in v.keys():
                if len(v[tag]) > int(settings["downsample_size"]):
                    v[tag] = set(random.sample(v[tag], int(settings["downsample_size"])))
    

        # Filter seqs to only contain the subsample
        keys_to_delete = [ ]
        for k in seqs.keys():
            new_entry = []
            for s in seqs[k]:
                if s[3][0] in tag_mapping \
                  and s[3][1] in tag_mapping[s[3][0]] \
                  and s[3][2] in tag_mapping[s[3][0]][s[3][1]]:
                  new_entry.append(s)
            seqs[k] = new_entry
            if not len(seqs[k]):
                keys_to_delete.append(k)
        for k in keys_to_delete:
            del seqs[k]
        print "Number of unique genotypes after downsampling: %d" % len(seqs)
        
        
    
    # Read reference sequences if given in config.
    if settings["reference_seqs"] != "none":
        print "Adding reference sequences to database..."
        ref_cnt = 0
        for r_f in settings["reference_seqs"].split():
            ref_fasta = SeqIO.parse(open(r_f), "fasta")

            for seq_record in ref_fasta:
                if len(seq_record.seq) == 0:
                    continue
                ref_cnt += 1
                seqs[str(seq_record.seq).upper()] = \
                  seqs.get(str(seq_record.seq).upper(), [] ) + [["_____ref", seq_record.id, False]]

        print "Number of reference sequences: ", ref_cnt


    print "Total number of genotypes in analysis: ", len(seqs.keys())

    if len(seqs.keys()) == 0:
	ut = open(settings["output_dir"] + "/no_sequences_to_cluster.txt", "wct")
	ut.write("Cowardly refusing to cluster 0 sequences!")
	ut.close()
	return


    # Collapse homopolymers if requested, otherwise just create data structure.

    new_seqs = dict()

    for seq in seqs.keys():
        if "max_homopolymer" in settings and int(settings["max_homopolymer"]) > 0:
            new_seq = reduce_homopolymer(seq,int(settings["max_homopolymer"]))
        else:
            new_seq = seq
            
        if new_seq not in new_seqs:
            new_seqs[new_seq] = dict()
        new_seqs[new_seq][seq]=seqs[seq]
        
    seqs = new_seqs
                                      
        
    print "Total number of unique seqs after homopolymer collapse: ", len(seqs.keys())

    uniseq_to_seq = UniseqDB(settings["work_dir"] + "/uniseq_to_seq.pick", "w")


    # Count sequences and build uniseq_to_seq

    for cs in enumerate(seqs):

        uniseq = {"seq"  : cs[1],
                  "tags" : dict(),
                  "revs" : 0,
                  "count" : 0,
                  "seqs" : dict()}

        for ucs in enumerate(seqs[cs[1]]):
            if ucs[1] not in uniseq["seqs"]:
                uniseq["seqs"][ucs[1]] = dict(count=0,
                                              revs=0,
                                              tags=dict())
                
            for s_info in seqs[cs[1]][ucs[1]]:
                if s_info[0] != "_____ref":
                    uniseq["count"] += 1
                    uniseq["seqs"][ucs[1]]["count"] += 1
                    if s_info[2]:
                        uniseq["revs"] += 1
                        uniseq["seqs"][ucs[1]]["revs"] += 1
                        
                if s_info[0] not in uniseq["tags"]:
                    uniseq["tags"][s_info[0]] = list()
                uniseq["tags"][s_info[0]].append(s_info)

                if s_info[0] not in uniseq["seqs"][ucs[1]]["tags"]:
                    uniseq["seqs"][ucs[1]]["tags"][s_info[0]] = list()
                uniseq["seqs"][ucs[1]]["tags"][s_info[0]].append(s_info)
                
        uniseq_to_seq["seq" + str(cs[0])] = uniseq


    # Close shelve to make sure it is written to disk
    uniseq_to_seq.sync()

    cPickle.dump(uniseq_to_seq.keys(),open(settings["work_dir"] + "/seqs.pick",
                                  "w"))


    # Prepare blast database


    # Run sub-clustering using som grid middleware (sge) for parallellisation
    if pr:
        pr.set_msg("(1/3) Clustering: ")
        pr.reset()

    cluster_program = "cluster_blast.py"
    mem_need="1G"

    step = 4000
    
    if "cluster_engine" in settings:
        if settings["cluster_engine"] == "usearch":
            cluster_program = "cluster_usearch.py"
            step = 4000
        elif settings["cluster_engine"] == "vsearch":
            cluster_program = "cluster_vsearch.py"
            step = 10000
	   
    

    
    num_jobs = len(seqs) / step + 1 if len(seqs) % step else len(seqs) / step

    files = [ ]

    job_num=0;


    print "Preparing Grid job"
    xgj_mem=1
    xgj = GridJob("ScataC" + config_file,
                  settings["work_dir"], settings["sge_params"],
                  59, str(xgj_mem) + "G", pr )


    for i in range(0,num_jobs):
        for j in range(i,num_jobs):
            job_num += 1
            xgj.add_task( (settings["lib_dir"] + "/" + cluster_program), 
                          [ settings["work_dir"], str(i), str(j), str(step), str(job_num) + ".fas" ])
            files.append(str(job_num) + ".fas")
    print "Number of blast jobs: ", len(files)

    xgj.start()

    print "Waiting for Grid clustering job to finish"
    xgj.check()
    xgj.wait()

    while xgj.get_num_failed() > 0:
        print "Some jobs failed, retrying"
        xgj.reset_failed()
        xgj_mem = xgj_mem * 2
        xgj.vf=str(xgj_mem) + "G"

        if pr:
            pr.set_msg("(1/3) Clustering (retry): ")
            pr.reset()
        xgj.start()
        xgj.check()
        xgj.wait()

    cpu += xgj.get_time()["cpu"]
    wallclock += xgj.get_time()["wallclock"]

    print "Clustering done"

    xgj.clean()


    # Merge clusters from all cross-blast runs.


    # This is done by hierarchical joining
    xgj2_mem=4
    retries = 0
    while True:
        retries += 1
        try:
            print "Merging clusters..."
            pick_files = map(lambda x: [settings["work_dir"] + "/" + x + ".pick", []], files)
            #print pick_files
            xgj2 = GridJob("ScataM" + config_file,
                           settings["work_dir"], settings["sge_params"], 90, str(xgj2_mem) + "G", pr)

            outnum = 1
            while True:
                partitioned_pick = [pick_files[x:x+50] 
                                    for x in range(0, len(pick_files), 50)]
                pick_files = [ ]

                for p in partitioned_pick:
                    filename = settings["work_dir"] + "/merge_out" + str(outnum) + ".pick"
                    id = xgj2.add_task(settings["lib_dir"] + "/merge_clusters.py",
                                       [ settings["work_dir"], filename ] +
                                       map(lambda x: x[0], p), 
                                       reduce(lambda a,b: a + b, map(lambda x: x[1], p)))
                    pick_files.append([filename, [id]])
                    outnum += 1

                if len(pick_files) == 1:
                    #print pick_files
                    break


            if pr:
                pr.set_msg("(2/3) Merging: ")
                pr.reset()

            xgj2.start()
            xgj2.wait()
            while xgj2.get_num_failed() > 0:
                print "Some jobs failed, retrying"
                xgj2.reset_failed()
                xgj2_mem = xgj2_mem * 2
                xgj2.vf=str(xgj2_mem) + "G"

                if pr:
                    pr.set_msg("(2/3) Retry merging: ")
                    pr.reset()
                xgj2.start()
                xgj2.check()
                xgj2.wait()

            cpu += xgj2.get_time()["cpu"]
            wallclock += xgj2.get_time()["wallclock"]


            clusters = map(lambda a: dict( set = a ), 
                           cPickle.load(open(pick_files[0][0])))
            if int(settings["graph"]) == 1:
                links = cPickle.load(open(pick_files[0][0] + ".links"))
            xgj2.clean()
            break
        except OSError:
            if retries > 7:
                raise Exception("Too many retries")
            print "Retrying"

    while True:
        try:

            # Go through clusters, assign names and if not singleton, add to 
            # job list for summarising cluster stats
            print "Summarising clusters"
            clusters_to_summarise = [ ]

            for c in enumerate(clusters):
                cl = c[1]
                num = c[0]

                cl["id"] = settings["job_id"] + "_" + str(num)
                cl["num"] = num

                if len(cl["set"]) > 1 or len(uniseq_to_seq[list(cl["set"])[0]]["seqs"]) > 1:
                    clusters_to_summarise.append(num)
                else:
                    refs = ""
                    if '_____ref' in uniseq_to_seq[list(cl["set"])[0]]["tags"]:
                        refs_list = uniseq_to_seq[list(cl["set"])[0]]["tags"]["_____ref"]
                        refs_list = [a[1] for a in refs_list]
                        refs = ", ".join(refs_list)

                    repseqs = [ ]
                    for uniseq in cl["set"]:
                        for subseq in uniseq_to_seq[uniseq]["seqs"].keys():
                            repseqs.append([uniseq, subseq, uniseq_to_seq[uniseq]["seqs"][subseq]["count"]])



                    repseqs.sort(lambda a, b: cmp(b[2], a[2]))
    
                    repseq = uniseq_to_seq[repseqs[0][0]]["seq"]
                    repseqs += [["","",0]]* int(settings["num_repseqs"]) # Padding
                    repseqs = repseqs[0:int(settings["num_repseqs"])]
                    
                    cl["cl_stats"] = dict( repseqs = repseqs,
                                           total_size = uniseq_to_seq[list(cl["set"])[0]]["count"],
                                           rev_count  = uniseq_to_seq[list(cl["set"])[0]]["revs"],
                                           proportion_singletons = 0,
                                           homop_factor = len(uniseq_to_seq[list(cl["set"])[0]]["seqs"]),
                                           genotypes = 1,
                                           references = refs,
                                           foo = True)


            cPickle.dump(clusters,open(settings["work_dir"] + "/clusters.pick",
                                      "w"))
            #cPickle.dump(uniseq_to_seq,open(settings["work_dir"] + "/seq_to_uniseq.pick",
            #                               "w"))
            xgj3_mem=3
            xgj3 = GridJob("ScataS" + config_file,
                           settings["work_dir"], settings["sge_params"],
                           30, str(xgj3_mem) + "G", pr )

            cts = [str(a) for a in clusters_to_summarise]
            for cl_group in [cts[x:x+100] \
                             for x in range(0, len(cts), 100)]:
                xgj3.add_task(settings["lib_dir"] + "/summarise_cluster.py",
                              [ settings["work_dir"]] + cl_group,[])

            xgj3.start()
            if pr:
                pr.set_msg("(3/3) Summarising: ")
                pr.reset()

            xgj3.check()
            xgj3.wait()

            while xgj3.get_num_failed() > 0:
                print "Some jobs failed, retrying"
                xgj3.reset_failed()
                xgj3_mem = xgj3_mem * 2
                xgj3.vf=str(xgj3_mem) + "G"
                xgj3.start()
                if pr:
                    pr.set_msg("(3/3) Retry summarising: ")
                    pr.reset()

                xgj3.check()
                xgj3.wait()

            cpu += xgj3.get_time()["cpu"]
            wallclock += xgj3.get_time()["wallclock"]

            for cluster in clusters_to_summarise:
                d = cPickle.load(open(settings["work_dir"] + "/cluster" + str(cluster) + ".pick"))
                clusters[cluster]["cl_stats"] = d

            #cPickle.dump(clusters,open(settings["work_dir"] + "/clusters.pick",
            #                              "w"))
            #if int(settings["graph"]) == 1:
            #    cPickle.dump(links,open(settings["work_dir"] + "/links.pick",
            #                           "w"))


            xgj3.clean()
            break
        except ValueError:
            print "Retrying"

    if pr:
    	pr.set_msg("Done.")
	pr.update_state(1,1)

    # Count singletons

    singleton_cnt = 0
    identified_clusters = 0
    identified_singletons = 0
    foo = 0
    foo_c = 0
    foo_s = 0
    bar = 0
    bar_c = 0
    bar_s = 0
    
    for c in clusters:
        if "foo" in c["cl_stats"]:

            foo += c["cl_stats"]["total_size"]
            if  c["cl_stats"]["total_size"] == 1:
                foo_s += 1
                #print c
            else:
                foo_c += 1
        else:
            bar += c["cl_stats"]["total_size"]
            if  c["cl_stats"]["total_size"] == 1:
                bar_s += 1
            else:
                bar_c += 1

        #print
        #print c["id"], c["cl_stats"]["total_size"]
        #for u in list(c["set"]):
            #print uniseq_to_seq[u]

        
        if c["cl_stats"]["total_size"] == 1:
            singleton_cnt += 1
            if len(c["cl_stats"]["references"]) > 0:
                identified_singletons += 1
        else:
            if len(c["cl_stats"]["references"]) > 0:
                identified_clusters += 1


    # Sort the clusters according to size

    clusters.sort(None,lambda cl: cl["cl_stats"]["total_size"],True)

    # Create a summary file and a fasta for all clusters
    print "Creating general cluster summary statistics"

    cluster_summary = open(settings["output_dir"] +
                           "/all_clusters_" + settings["job_id"] + ".txt","w")
    cluster_singleton_summary = open(settings["output_dir"] +
                                     "/all_singletons_" + settings["job_id"] +
                                     ".txt","w")
    cluster_fas = open(settings["output_dir"] + "/all_clusters_" +
                       settings["job_id"] + ".fas","w")
    cluster_singleton_fas = open(settings["output_dir"] +
                                 "/all_singletons_" + settings["job_id"] + ".fas","w")

    cluster_summary.write("General statistics for clustering run " +
                          settings["job_name"] + "\n")



    to_print=[ ["pyro_reads", "Total number of pyro reads"],
               ["good_reads",  "Sequences passing QC"],
               ["qual_type", "Quality screening method"],               
               ["rev", "Reads matching after reverse complement"],
               ["truncated", "Number of reads truncated"],
               ["total_skipped", "Number of reads discarded"],
               ["too_short", "Reads too short"],
               ["low_mean_quality", "Reads with too low mean quality"],
               ["low_min_quality", "Reads containing bases with too low quality"],
               ["mean_len", "Mean readlength"],
               ["max_len", "Maximum read length"],
               ["min_len", "Minimum read length"],
               ["no_primer5", "Missing 5' primer"],
               ["no_primer3", "Missing 3' primer"],
               ["no_tag5", "Missing 5' tag"],
               ["no_tag3", "Missing 3' tag"]]
    no_sum = set([ "max_len", "min_len", "mean_len", "qual_type" ])

    names = settings["454names"].split()

    if len(names) > 1:
        cluster_summary.write(";" + ";".join(names) + ";Total\n")
    else:
        cluster_summary.write(";" + ";".join(names) + "\n")

    for p in to_print:
        cluster_summary.write(p[1] + ":")
        s = 0
        for n in names:
            if p[0] not in no_sum:
                s += seq_stats[n][p[0]]
            if p[0] in seq_stats[n]:
                cluster_summary.write(";" + str(seq_stats[n][p[0]]))
            else:
                cluster_summary.write(";")
        if len(names) > 1:
            if p[0] not in no_sum:
                cluster_summary.write(";" + str(s) + "\n")
            else:
                cluster_summary.write(";\n")
        else:
            cluster_summary.write("\n")


    
    cluster_summary.write("""

Number of global clusters:;%d
Number of identified cluster:;%d
Number of global singletons:;%d
Number of identified singletons:;%d



Comparison parameters:
Clustering engine:; %s
Minimum alignment length for clustering:;%s
Maximum distance:;%f
Missmatch penalty:;%f
Gap open penalty:;%f
Gap extension penalty:;%f
End gap weight:;%f
Homopolymer reduction at:;%d

Resource usage:
Wall clock time:;%.3f
Core time:;%.3f


""" % ( 
            len(clusters) - singleton_cnt,
            identified_clusters,
            singleton_cnt,
            identified_singletons,
            settings["cluster_engine"],
            settings["min_alignment"],
            float(settings["max_dist"]),
            float(settings["missmatch_pen"]),
            float(settings["gap_open_pen"]),
            float(settings["gap_extend_pen"]),
            float(settings["end_gap_weight"]),
            int(settings["max_homopolymer"]),
            wallclock,
            cpu))

    cluster_summary.write("Cluster ID;Cluster Size;Num reversed;Reference;Proportion Singletons;Num genotypes;Homopolymer factor;")
    for i in range(0,int(settings["num_repseqs"])):
        cluster_summary.write("Genotype%d Freq;Sequence%d;" % (i + 1, i + 1))
    cluster_summary.write("\n")
    cluster_singleton_summary.write("Cluster ID;Cluster Size;Reference;Repseq\n")
    for c in clusters:
        #print c["cl_stats"]
        if c["cl_stats"]["total_size"] > 1:
            cluster_summary.write("%s;%d;%d;%s;%f;%d;%f;" % (
                    c["id"],
                    c["cl_stats"]["total_size"],
                    c["cl_stats"]["rev_count"],
                    c["cl_stats"]["references"],
                    c["cl_stats"]["proportion_singletons"],
                    c["cl_stats"]["genotypes"],
                    c["cl_stats"]["homop_factor"]
                ))
            for i in range(0, int(settings["num_repseqs"])):
                if c["cl_stats"]["repseqs"][i][2] > 0:
                    cluster_summary.write("%d;%s;" % (c["cl_stats"]["repseqs"][i][2], c["cl_stats"]["repseqs"][i][1]))
                else:
                    cluster_summary.write("0;;")
            cluster_fas.write(">" + c["id"] + "_repseq_" + str(0) + "\n" + c["cl_stats"]["repseqs"][0][1] + "\n")
            cluster_summary.write("\n")

        else:
            cluster_singleton_summary.write("%s;%d;%s;%f;%s\n" % (
                c["id"],
                c["cl_stats"]["total_size"],
                c["cl_stats"]["references"],
                c["cl_stats"]["proportion_singletons"],
                c["cl_stats"]["repseqs"][0][1]))
            cluster_singleton_fas.write(">" + c["id"] + "\n" + c["cl_stats"]["repseqs"][0][1] + "\n")


    cluster_summary.close()
    cluster_fas.close()

    print "Dumping cluster alignments"
    for c in filter(lambda c: c["cl_stats"]["total_size"] > 1, 
                    clusters):
        ali_fas = open(settings["output_dir"] + "/clusters/" + c["id"] + ".fas","w")
        for uniseq in list(c["set"]):
            for subseq in uniseq_to_seq[uniseq]["seqs"]:
                for (tag, ids) in uniseq_to_seq[uniseq]["seqs"][subseq]["tags"].iteritems():
                    #print tag, ids
                    for sid in ids:
                        ali_fas.write(">" + sid[0] + "_" + sid[1] + \
                                      ("_rev_comp" if sid[2] else "" ) + "\n" + subseq + "\n")

        if "repseq_alist" in c["cl_stats"]:
            for rseq in c["cl_stats"]["repseq_alist"]:
                ali_fas.write(">" + rseq["id"] + "\n" + rseq["seq"].replace("-", "") + "\n")

        ali_fas.close()
    print "Dumping cluster reference alignments"
    for c in filter(lambda c: c["cl_stats"]["total_size"] > 1, 
                    clusters):
	if "repseq_alist" in c["cl_stats"]:
	        ali_fas = open(settings["output_dir"] + "/repseq_alignments/" + c["id"] + ".fas","w")
                for rseq in c["cl_stats"]["repseq_alist"]:
                    ali_fas.write(">" + rseq["id"] + "\n" + rseq["seq"] + "\n")
        	ali_fas.close()

    print "Summarising results per tag"
    tag_stats = { }
    global_tag_summary = open(settings["output_dir"] + "/all_tags.txt","w")
    tag_by_cluster = open(settings["output_dir"] + "/all_tag_by_cluster.txt", "w")
    tag_by_cluster_c = open(settings["output_dir"] + "/all_tag_by_cluster_counts.txt", "w")


    global_tag_summary.write("Tag name;Size;Prop of global;Singletons;Prop singletons;Global singletons;Prop global;Prop of total global;#1;#2;#3;#4;#5;#6;#7;#8;#9;#10\n")

    tag_by_cluster_clusters = filter(lambda c: c["cl_stats"]["total_size"] > 1, 
                                     clusters)[:int(settings["tag_by_cluster_max"])]

    tag_by_cluster.write("Tag")
    tag_by_cluster_c.write("Tag")
    for c in tag_by_cluster_clusters:
        tag_by_cluster.write(";" + c["id"])
        tag_by_cluster_c.write(";" + c["id"])
    tag_by_cluster.write("\n")
    tag_by_cluster_c.write("\n")


    for tag in used_tags.iterkeys():

        #print "Doing tag %s..." % (tag)
        tag_stats[tag] = dict(total_size = 0,
                              global_singletons = 0,
                              clusters = list(),
                              clusters_d = dict())

        for i, c in enumerate(clusters):
            tmp_stat = dict()
            tmp_stat["c_ind"] = i        
            tmp_stat["cnt"] = 0
            tmp_stat["references"] = c["cl_stats"]["references"]
            for seq in c["set"]:
                if tag in uniseq_to_seq[seq]["tags"]:
                    tag_stats[tag]["total_size"] += \
                        len(uniseq_to_seq[seq]["tags"][tag])
                    tmp_stat["cnt"] += len(uniseq_to_seq[seq]["tags"][tag])
            if tmp_stat["cnt"] == 1 and c["cl_stats"]["total_size"] == 1:
                tag_stats[tag]["global_singletons"] += 1
                tmp_stat["global_singleton"] = 1
            else:
                tmp_stat["global_singleton"] = 0
            if tmp_stat["cnt"] >= 1:
                tmp_stat["seq"] = c["cl_stats"]["repseqs"][0][1]
                tmp_stat["id"] = c["id"]
                tag_stats[tag]["clusters"].append(tmp_stat)
        tag_stats[tag]["singletons"] = filter(lambda a: a["global_singleton"] == 1,
                                              tag_stats[tag]["clusters"])
        tag_stats[tag]["clusters"] = filter(lambda a: a["global_singleton"] != 1,
                                            tag_stats[tag]["clusters"])
        # Sort cluster descending
        tag_stats[tag]["clusters"].sort(None, lambda a: a["cnt"], True)


        # Make cluster dict for later tag_by_cluster output
        for c in tag_stats[tag]["clusters"]:
            tag_stats[tag]["clusters_d"][c["id"]] = c

        # Save stats for this tag into a file
        tag_summary = open(settings["output_dir"] + "/tag_summaries/" + tag.replace("/","_").replace(":","_") + "_clusters.txt","w")
        tag_singleton_summary = open(settings["output_dir"] + "/tag_summaries/" + tag.replace("/","_").replace(":","_") + "_singletons.txt","w")
        tag_fas = open(settings["output_dir"] + "/tag_fastafiles/" + tag.replace("/","_").replace(":","_") + "_clusters.fas","w")
        tag_singleton_fas = open(settings["output_dir"] + "/tag_fastafiles/" + tag.replace("/","_").replace(":","_") + "_singletons.fas","w")            

        tag_summary.write("""Summary for tag %s

    Total number of sequences:; %d
    Proportion of global sequences:; %f
    Number of singletons:; %d
    Number of global singletons:; %d
    Proportion of singletons being global singletons:; %f
    Proportion of global singletons:; %f


    Clusters:
    Cluster id;Reference;Cluster size:Repseq
    """ %  (tag, 
            tag_stats[tag]["total_size"],
            tag_stats[tag]["total_size"] / float(global_seqs),
            len(tag_stats[tag]["singletons"]),
            tag_stats[tag]["global_singletons"],
            tag_stats[tag]["global_singletons"] / \
                float(len(tag_stats[tag]["singletons"])) if len(tag_stats[tag]["singletons"]) else -1.0,
            (tag_stats[tag]["global_singletons"] / float(singleton_cnt)) if singleton_cnt else -1.0,
            ))
        global_tag_summary.write("%s;%d;%f;%d;%f;%d;%f;%f" % 
                                 (tag, 
                                  tag_stats[tag]["total_size"],
                                  tag_stats[tag]["total_size"] / \
                                      float(global_seqs),
                                  len(tag_stats[tag]["singletons"]),
                                  (len(tag_stats[tag]["singletons"]) / \
                                      float(tag_stats[tag]["total_size"])) if tag_stats[tag]["total_size"] else -1.0,
                                  tag_stats[tag]["global_singletons"],
                                  tag_stats[tag]["global_singletons"] / \
                                      float(len(tag_stats[tag]["singletons"])) if len(tag_stats[tag]["singletons"]) \
                                                                                  else -1.0,
                                  tag_stats[tag]["global_singletons"] / float(singleton_cnt) if singleton_cnt else -1.0,
                                  ))

        tag_singleton_summary.write("Cluster id;Reference;Cluster size:Repseq\n")
        cnt = 0
        for c in tag_stats[tag]["clusters"]:
            tag_summary.write("%s;%s;%d;%s\n" % (c["id"], c["references"], 
                                                 c["cnt"], c["seq"]))
            tag_fas.write(">%s\n%s\n" % (c["id"], c["seq"]))
            if cnt < 10:
                global_tag_summary.write(";%d" % (c["cnt"]))
            cnt += 1
        global_tag_summary.write("\n")
        for c in tag_stats[tag]["singletons"]:
            tag_singleton_summary.write("%s;%s;%d;%s\n" % (c["id"], c["references"], c["cnt"], c["seq"]))
            tag_singleton_fas.write(">%s\n%s\n" % (c["id"], c["seq"]))


        tag_by_cluster.write(tag)
        tag_by_cluster_c.write(tag)
        for c in tag_by_cluster_clusters:
            if c["id"] in tag_stats[tag]["clusters_d"]:
                t=tag_stats[tag]["clusters_d"][c["id"]]
                if tag_stats[tag]["total_size"]:
                    tag_by_cluster.write(";" + str(float(t["cnt"]) / float(tag_stats[tag]["total_size"])))
                    tag_by_cluster_c.write(";" + str(float(t["cnt"])))
                else:
                    tag_by_cluster.write(";NA")
                    tag_by_cluster_c.write(";NA")
            else:
                tag_by_cluster.write(";0")
                tag_by_cluster_c.write(";0")
        tag_by_cluster.write("\n")
        tag_by_cluster_c.write("\n")
        #break

    #print tag_stats

    if int(settings["graph"]) == 1:
        print "Creating gv graph files for all clusters."
        for c in clusters:
            #print c
            cl_out = open(settings["output_dir"] + "/" + c["id"] + ".gv", "w")
            if c["cl_stats"]["total_size"] > 1:
                #print c["id"] + "..."
                cl_out.write("graph G {\n")

                for s in list(c["set"]):
                    ref = ""

                    if "_____ref" in uniseq_to_seq[s]["tags"]:
                        ref = "\\n" +", ".join(uniseq_to_seq[s]["tags"]["_____ref"])
                    cl_out.write('%s [ label = "%s\\n%d%s" ];\n' % 
                                 (s, s, uniseq_to_seq[s]["count"],ref))
                my_links = filter(lambda l: l["s1"] in c["set"] or 
                                  l["s2"] in c["set"], links)
                for link in my_links:
                    if link["mode"] == "current":
                        color = "black"
                    elif link["mode"] == "new":
                        color = "green"
                    elif link["mode"] == "join":
                        color = "red"
                    else:
                        color = "blue"
                    cl_out.write("%s -- %s [ color = %s ];\n" % 
                                 (link["s1"], link["s2"], color))
                cl_out.write(
                    '''overlap = False;
    splines = True;
    label = "%s: clustering distance %s, total size %d"
    }
    ''' 
                    % (c["id"], settings["max_dist"], c["cl_stats"]["total_size"]))
                                                                                    #                                           c["cl_stats"]
    #             cluster_summary.write("%s;%d;%s;%f;%f;%f;%s\n" % (
    #                     c["id"],
    #                 c["cl_stats"]["total_size"],
    #                 ",".join(c["cl_stats"]["references"]),
    #                 c["cl_stats"]["mean_diff"],
    #                 c["cl_stats"]["max_diff"],
    #                 c["cl_stats"]["proportion_singletons"],
    #                 c["cl_stats"]["consensus"]))



