#!/usr/bin/env python

import constants
import sys, time, sge, re, os, traceback, difflib, cPickle, subprocess
from Bio.Seq import Seq
from Bio.Alphabet import generic_dna
from array import array



sys.path.append("/mykopat/Linux-x86_64/lib/python2.5/site-packages/")
import MySQLdb

import PyroParser, FastqParser, FilterSeq

from constants import *


def master_loop(argv, mpi_checker):
    datasetid = int(argv[1])

    log_entry("Starting %s %s" % (argv[0], argv[1]))

    log_entry("Connecting to database.")

    print argv

    db = MySQLdb.connect( host=db_host, user=db_user, passwd=db_pass, db=db_db)
    db_c = db.cursor(MySQLdb.cursors.DictCursor)


    db_c.execute("""SELECT d.*, IF(ISNULL(t5.tagSetName),"(none)",t5.tagSetName) AS t5n,
                                IF(ISNULL(t3.tagSetName),"(none)",t3.tagSetName) t3n
                                FROM Datasets d
                                LEFT JOIN TagSets as t5 ON (d.Tagset5 = t5.idTagSets)
                                LEFT JOIN TagSets as t3 on (d.Tagset3 = t3.idTagSets)
                            WHERE idDatasets = %s""", (datasetid,))

    row = db_c.fetchone()

    if not row:
        if mpi_checker:
                mpi_checker.end_workers()
        log_entry("Job is gone!")
        sys.exit(0)

    owner = row["owner"]
    name = row["Name"]
    file_type=row["file_type"]

    if row["ready"] == 1:
        log_entry("Dataset %d already set to ready" % (datasetid))
        return 0
    

    db_c.execute("SELECT Email FROM Users WHERE idUsers = %s", (owner,))
    email = db_c.fetchone()["Email"]

    log_entry("Checking dataset %d for user %s." \
              % (datasetid, email))

    tags5 = get_tagset(row["Tagset5"])
    if len(tags5) == 0:
        print "No 5' tags"
    tags3 = get_tagset(row["Tagset3"])
    if len(tags3) == 0:
        print "No 3' tags"



    lengths = [ ]
    seqs = dict()
    msg = ""
    errors = ""
    try:
        if file_type == "fasta":
            qual_seqs = PyroParser.RawPyroRes(dataset_dir + "/" + sys.argv[1] + ".1.dat")
        elif file_type == "qual":
            qual_seqs = PyroParser.RawPyroRes(dataset_dir + "/" + sys.argv[1] + ".1.dat",
                                              dataset_dir + "/" + sys.argv[1] + ".2.dat")
        elif file_type == "sff":
            # HQR implies no clipping, full sequence implies clipping
            extract_status = subprocess.call([ base_dir + "/sff_extract",
                            "-u" if row["raw_filtering"] else "-c",
                            "-s", dataset_dir + "/" + sys.argv[1] + ".fasta",
                            "-q", dataset_dir + "/" + sys.argv[1] + ".qual",
                            dataset_dir + "/" + sys.argv[1] + ".1.dat" ])
            if extract_status:
                raise Exception("Unable to read sff file!")
            qual_seqs = PyroParser.RawPyroRes(dataset_dir + "/" + sys.argv[1] + ".1.fasta",
                                              dataset_dir + "/" + sys.argv[1] + ".2.qual")
        elif file_type == "fastq":
            qual_seqs = FastqParser.Single(dataset_dir + "/" + sys.argv[1] + ".1.dat")
        elif file_type == "fastq2":
            qual_seqs = FastqParser.Pair(dataset_dir + "/" + sys.argv[1] + ".1.dat",
                                         dataset_dir + "/" + sys.argv[1] + ".2.dat",
                                         row["overlap_kmer"], row["overlap_hsp"],
                                         row["overlap_min"])

        seq_stats = dict( pyro_reads = 0,
                          good_reads = 0,
                          unique_reads = 0,
                          total_skipped = 0,
                          too_short = 0,
                          truncated = 0,
                          low_mean_quality = 0,
                          low_min_quality = 0,
                          shorter_than_primer = 0,
                          no_primer5 = 0,
                          no_primer3 = 0,
                          no_tag5 = 0,
                          no_tag3 = 0,
                          rev = 0)

        # Grab 100 reads at a time
        processed_reads = [ ]
        chunk_size = 200
        mpi_listners = []

        while True:
            reads = [ ]
            print "grabbing more"
            for rread in qual_seqs:
                reads.append(rread)
                if len(reads) >= chunk_size:
                    break
            read_job = dict( reads = reads,
                             mean_qual = row["mean_qual"],
                             min_len = row["min_len"],
                             min_qual = row["min_qual"],
                             max_len = row["max_len"],
                             raw_filtering = row["raw_filtering"],
                             p5 = row["Primer5"],
                             p5s = row["Primer5score"],
                             p3 = str(Seq(row["Primer3"],generic_dna).reverse_complement()),
                             p3s = row["Primer3score"],
                             tags5 = tags5,
                             tags3 = tags3)
            if mpi_checker:
                while mpi_checker.data_available():
                    t = mpi_checker.get_data()
                    processed_reads += t
                    print len(processed_reads)
            if mpi_checker and mpi_checker.ready_to_send():
                mpi_checker.send(read_job)
            else:
                print "Processing in master"
                processed_reads += process_reads(read_job)
                print len(processed_reads)
            if len(reads) < chunk_size:
                break
        if mpi_checker:
            while not mpi_checker.no_more():
                while mpi_checker.data_available():
                    processed_reads += mpi_checker.get_data()
                    print len(processed_reads)
                time.sleep(1)
            mpi_checker.end_workers()
                
        print len(processed_reads), processed_reads[0:10]
        for r in processed_reads:
            seq_stats["pyro_reads"] += 1
            for s in r["status"]:

                if s in seq_stats:
                    seq_stats[s] += 1
                else:
                    seq_stats[s] = 1
            if "keep" in r:
                seq_stats["good_reads"] += 1
                detagged_seq = r["detagged_seq"]
                seq_name = detagged_seq["tag_name"] + "_" + r["id"]

                try:
                    seqs[detagged_seq["seq"]].append([(row["Name"] + " " + \
                                                   detagged_seq["tag_name"]),
                                                   seq_name, detagged_seq["rev"],
                        [row["Name"], detagged_seq["tag_name"], r["id"]]])
                except KeyError:
                    seqs[detagged_seq["seq"]] = [[(row["Name"] + " " + \
                                               detagged_seq["tag_name"]),
                                               seq_name,
                                               detagged_seq["rev"],
                        [row["Name"], detagged_seq["tag_name"], r["id"]]]]

                if len(detagged_seq["seq"]) == 0:
                    print detagged_seq

                lengths.append(len(detagged_seq["seq"]))
                if detagged_seq["rev"]:
                    seq_stats["rev"] += 1

        seq_stats["total_skipped"] = seq_stats["pyro_reads"] - seq_stats["good_reads"]
        seq_stats["unique_reads"] = len(seqs)
        seq_stats["mean_len"] = (sum(lengths) / seq_stats["good_reads"]) if seq_stats["good_reads"] > 0 else -1
        seq_stats["max_len"] = max(lengths) if (len(lengths) > 0) else -1
        seq_stats["min_len"] = min(lengths) if (len(lengths) > 0) else -1
        seq_stats["tagset5name"] = row["t5n"]
        seq_stats["tagset3name"] = row["t3n"]
        if row["raw_filtering"] == 0:
            seq_stats["qual_type"] = "Full sequence"
            if not qual_seqs.qual_present:
                seq_stats["qual_type"] += " no quality data provided"
        elif row["raw_filtering"] == 1:
            seq_stats["qual_type"] = "Full sequence - quality data ignored"
        elif row["raw_filtering"] == 2:
            seq_stats["qual_type"] = "HQR"
        elif row["raw_filtering"] == 3:
            seq_stats["qual_type"] = "Amplicon quality"
        else:
            seq_stats["qual_type"] = "Unknown filtering?"


        lengths.sort()
        print seq_stats
        cPickle.dump(seqs,
                     open(dataset_dir + "/" + sys.argv[1] + ".seqs.pick", "wct"))
        cPickle.dump(seq_stats,
                     open(dataset_dir + "/" + sys.argv[1] + ".stat.pick", "wct"))

        to_print=[ ["pyro_reads", "Total number of reads"],
                   ["good_reads",  "Sequences passing QC"],
                   ["qual_type", "Quality screening method"],
                   ["rev", "Reads matching after reverse complement"],
                   ["truncated", "Number of reads truncated (limit %d)" % (row["max_len"]) ],
                   ["total_skipped", "Number of reads discarded"],
                   ["failed_pair", "Failed overlap pairing (if applicable)"],
                   ["too_short", "Reads too short (< %d bp)" % ( row["min_len"] )],
                   ["low_mean_quality", "Reads with too low mean quality ( < %d)" % (row["mean_qual"])],
                   ["low_min_quality", "Reads containing bases with too low quality ( < %d )" % (row["min_qual"])],
                   ["mean_len", "Mean readlength"],
                   ["max_len", "Maximum read length"],
                   ["min_len", "Minimum read length"],
                   ["no_primer5", "Missing 5' primer (%s at %f)" % (row["Primer5"], row["Primer5score"]) ],
                   ["no_primer3", "Missing 3' primer (%s at %f)" % (row["Primer3"], row["Primer3score"])],
                   ["tagset5name", "Tagset 5"],
                   ["no_tag5", "Missing 5' tag"],
                   ["tagset3name", "Tagset 3"],
                   ["no_tag3", "Missing 3' tag"],
                   ["shorter_than_primer", "High quality sequence shorter than primer"] ]

        for tp in to_print:
            msg += "%s:%s%8s\n" % (tp[1], (" " * (54 - len(tp[1]))), str(seq_stats[tp[0]]) if tp[0] in seq_stats else "0")

    
        print msg
    
        if seq_stats["good_reads"] == 0:
            errors = "No good reads!\n\n" + msg 


    except:
	if mpi_checker:
		mpi_checker.end_workers()
        exceptionType, exceptionValue, exceptionTraceback = sys.exc_info()
        exception_text = "".join(traceback.format_exception(exceptionType, exceptionValue,
                                    exceptionTraceback))
        open(log_dir + "/check_dataset_exceptions_%04d.txt" % (datasetid),
             "w").write("%s\n%s" % (repr(sys.exc_info()), exception_text))
        send_mail("mikael.durling@mykopat.slu.se",
                  "Check DataSet %d failed" % (datasetid), 
                  "Owner %s\nDatasetid: %d\nException: %s\nBacktrace %s" \
                  % (email,
                     datasetid,
                     repr(sys.exc_info()),
                      exception_text))
        errors += "Parse error, please make sure the fasta and quality files are properly formated and in sync?\n"
        errors += "The parser reported: %s\n" % (exception_text)
        if mpi_checker:
            mpi_checker.end_workers()

    subject = "SCATA: dataset " + name
    body = "The dataset you submitted to SCATA has been processed.\n"
    if errors != "":
        log_entry("Errors found in dataset %d for user %s" % (datasetid, email))
        subject += " contains errors"
        body += "Errors were found while processing the dataset. An attempt to give details of the errors is found below. \n\n"
        body += "The dataset has been removed from the database. Pleas correct any errors and retry your upload.\n\n"
        body += "Error log:\n" + errors

        # Make sure we have a recent connection to the database.
        db = MySQLdb.connect( host=db_host, user=db_user, passwd=db_pass, db=db_db)
        db_c = db.cursor(MySQLdb.cursors.DictCursor)
        #TODO db_c.execute("DELETE FROM Datasets WHERE idDatasets = %s", (datasetid,))
        db.commit()
        try:
            print "not removing files"
            #os.remove(dataset_dir + "/" + sys.argv[1] + ".dat.1")
            #os.remove(dataset_dir + "/" + sys.argv[1] + ".dat.2")
        except:
            log_entry("Could not remove file.")
    else:
        log_entry("Dataset %d accepted for user %s. [%d sequences, mean length %d, max length %d]" \
                  % (datasetid, email, seq_stats["good_reads"],
                     seq_stats["mean_len"], seq_stats["max_len"]))
        subject += " is accepted"
        body += "No errors were found and the dataset can now be used in subsequent analyses.\n\n"
        body += "Below follows a summary of the dataset.\n\n%s\n\n" % \
            (msg)
        db = MySQLdb.connect( host=db_host, user=db_user, passwd=db_pass, db=db_db)
        db_c = db.cursor(MySQLdb.cursors.DictCursor)
        #db_c.execute("UPDATE Datasets SET ready=1,Description=%s WHERE idDatasets = %s", (msg,datasetid,))
        db_c.execute("UPDATE Datasets SET ready=0,locked=0,Description=%s WHERE idDatasets = %s", (msg,datasetid,))
        db.commit()
        try:
            print "removeing files"
            os.remove(dataset_dir + "/" + sys.argv[1] + ".fas")
            os.remove(dataset_dir + "/" + sys.argv[1] + ".qual")
        except:
            log_entry("Could not remove file.")

    send_mail(email,subject,body)




def get_tagset(tagset):
    if int(tagset) > 0:
        tags = dict()
        tag_file = open( tagset_dir + "/" + str(tagset) + ".txt")
        tag_length = 0
        # Parse tag-file and perform sanity checks
        for line in tag_file:
            line_list = list()
            if line == "":
                continue
            if line[-1] == "\n":
                line = line[:-1]
            if line!="" and line[-1] == "\r":
                line = line[:-1]

            line_list = line.split(";")
            #print line_list
            if len(line_list) != 2:
                continue
            tags[line_list[1]] = line_list[0]
            if tag_length == 0:
                tag_length = len(line_list[1])

            if tag_length != len(line_list[1]):
                raise Exception("Different length on tags")
            
            tags["_____length"] = tag_length
        return tags
    return { }

# Function to process a batch of reads
def process_reads(read_job):

    result_list = [ ]
    dtseq = FilterSeq.DeTagSeq(read_job["p5"], read_job["p5s"],
                              read_job["p3"], read_job["p3s"],
                              read_job["tags5"], read_job["tags3"])
    for qseq in read_job["reads"]:
        r=qseq.get()
        res = dict(status = dict())
        seq = r[0]
        qual = r[1]
        detagged_seq = None
        if seq:
            res["id"] = seq.id
            if read_job["raw_filtering"] == 0: # Filter full
                if not qual:
                    raise Exception("Full sequence filtering requires quality data")
                
                seq, status = FilterSeq.filter_full(seq, qual, read_job["min_len"],
                                                    read_job["mean_qual"],
                                                    read_job["min_qual"])
            elif read_job["raw_filtering"] == 1: # Full sequence, no filtering
                status = None
            elif read_job["raw_filtering"] == 2: # HQR
                if not qual:
                    raise Exception("HQR filtering requires quality data")
                seq, status = FilterSeq.filter_hqr(seq, qual, read_job["min_len"],
                                                    read_job["mean_qual"],
                                                    read_job["min_qual"])
            elif read_job["raw_filtering"] == 3: # Primers first
                if not qual:
                    raise Exception("Amplicon quality filtering requires quality data")
                detagged_seq, qual = dtseq.detag_seq(seq,qual)
                if "status" in detagged_seq:
                    res["status"][detagged_seq["status"]] = 1
                    result_list.append(res)
                    continue

                seq, status = FilterSeq.filter_full(detagged_seq["seq_record"], qual, read_job["min_len"],
                                                    read_job["mean_qual"],
                                                    read_job["min_qual"])
            else:
                raise Exception("Unknown filter type")
        else:
            status = "failed_pair"

        if status:
            res["status"][status] = 1
            result_list.append(res)
            continue

        if len(seq.seq) < read_job["min_len"]:
            res["status"]["too_short"] = 1
            result_list.append(res)
            continue
        if not detagged_seq:
            detagged_seq = dtseq.detag_seq(seq)
        #detagged_seq = FilterSeq.old_detag_seq(seq, read_job["p5"], read_job["p5s"],
        #                      read_job["p3"], read_job["p3s"],
        #                      read_job["tags5"], read_job["tags3"]) 
        if "status" in detagged_seq:
            res["status"][detagged_seq["status"]] = 1
            result_list.append(res)
            continue
        
        if len(detagged_seq["seq"]) < read_job["min_len"]:
            res["status"]["too_short"] = 1
            result_list.append(res)
            continue
        
        if read_job["max_len"] > 0 and len(detagged_seq["seq"]) > read_job["max_len"]:
            res["status"]["truncated"] = 1
            detagged_seq["seq"] = detagged_seq["seq"][:read_job["max_len"]]
        if read_job["max_len"] < 0:
            res["status"]["truncated"] = 1
            detagged_seq["seq"] = detagged_seq["seq"][:read_job["max_len"]]

        res["detagged_seq"] = detagged_seq
        res["keep"] = True
        res["id"] = seq.id
        result_list.append(res)
    return result_list


class Mpi_Checker:
    def __init__(self, mpi_comm, buf_size=5000):
        self.buf_size=buf_size
        self.mpi_comm = mpi_comm
        self.size = mpi_comm.Get_size()
        self.reqs = {}

        for i in range(1,self.size):
            self.reqs[i] = dict(available=True,
                                  req=None,
                                  status=None,
                                  ready=False,
                                  buf=array('c', '\0' * (self.buf_size + 100)))

    # Checks if any worker is done and returns True if at least one is so.
    def data_available(self):
        for k in self.reqs.keys():
            if self.reqs[k]["req"]:
                if self.reqs[k]["req"].Test(self.reqs[k]["status"]):
                    self.reqs[k]["ready"] = True
                    print "Ready", k
                    return True
        return False

    def get_data(self):
        for k in self.reqs.keys():
            if self.reqs[k]["ready"]:
                print "r", k
                pickled = self.reqs[k]["buf"].tostring().split("\0")[0]
                while True and pickled[0]=="c":
                    buf=array('c', '\0' * (self.buf_size + 100))
                    self.mpi_comm.Recv(buf,source=k,tag=42)
                    tmp = buf.tostring().split("\0")[0]
                    if tmp[0]=='l':
                        pickled += tmp[1:]
                        break
                    pickled += tmp[1:]
                self.reqs[k]["available"]=True
                self.reqs[k]["ready"]=False
                    
                return cPickle.loads(pickled[1:])
        return None

    def send(self, job):
        for k in self.reqs.keys():
            if self.reqs[k]["available"]:
                pickled = cPickle.dumps(job)
                print "Send", k, "pickled len", len(pickled)
                splitted = ["c" + pickled[x:x + (self.buf_size-1)] for x in range(0, len(pickled), self.buf_size - 1)]
                splitted[-1] = 'l' + splitted[-1][1:]
                for s in splitted:
                    self.mpi_comm.Send(array('c',s), dest=k, tag=42)
                self.reqs[k]["available"] = False
                self.reqs[k]["buf"]=array('c', '\0' * (self.buf_size + 100))
                self.reqs[k]["req"] = self.mpi_comm.Irecv(self.reqs[k]["buf"], source=k, tag=42)
                self.reqs[k]["status"] = MPI.Status()
                return
                
    def ready_to_send(self):
        for k in self.reqs.keys():
            if self.reqs[k]["available"]:
                return True
        return False

    def no_more(self):
        for k in self.reqs.keys():
            if not self.reqs[k]["available"]:
                return False
        return True

    def end_workers(self):
        for k in self.reqs.keys():
            self.mpi_comm.Send(array('c','l' + cPickle.dumps(None)), dest=k, tag=42)

            
class Mpi_Worker:
    def __init__(self, mpi_comm,buf_size=5000):
        self.buf_size=buf_size
        self.mpi_comm = mpi_comm
        self.rank = mpi_comm.Get_rank()

    # Worker, will exit once it gets the close message
    def worker(self):
        print "Worker in", self.rank, "now running"
        while True:
            buf = array('c', '\0' * self.buf_size)
            self.mpi_comm.Recv(buf, source=0, tag=42)
            pickled = buf.tostring().split('\0')[0]
            if pickled[0] == 'c':
                while True:
                    buf = array('c', '\0' * self.buf_size)
                    self.mpi_comm.Recv(buf, source=0, tag=42)
                    tmp = buf.tostring().split("\0")[0]
                    if tmp[0]=='l':
                        pickled += tmp[1:]
                        break
                    pickled += tmp[1:]
            print "Recv", self.rank, "pickled len", len(pickled)
            job = cPickle.loads(pickled[1:])
            if job == None:
                return
            result = process_reads(job)
            pickled = cPickle.dumps(result)
            splitted = ["c" + pickled[x:x + (self.buf_size - 1)] for x in range(0, len(pickled), self.buf_size - 1)]
            splitted[-1] = 'l' + splitted[-1][1:]
            for s in splitted:
                self.mpi_comm.Send(array('c',s), dest=0, tag=42)
        print "Worker in", self.rank, "now done"


mpi_checker = None
try:
    from mpi4py import MPI
    mpi_comm = MPI.COMM_WORLD
    mpi_rank = mpi_comm.Get_rank()
    mpi_size = mpi_comm.Get_size()
    print "MPI Available"
    if mpi_size == 1:
        print "Only one MPI process, will not use MPI"
        master_loop(sys.argv, None)
    else:
        if mpi_rank > 0:
            print "Starting MPI subprocess listener #", mpi_rank
            Mpi_Worker(mpi_comm).worker()
            print "MPI processor done #", mpi_rank
        else:
            print "MPI main process with ", mpi_size, "subprocesses"
            master_loop(sys.argv, Mpi_Checker(mpi_comm))
except ImportError:
    print "No MPI"
    master_loop(sys.argv, None)

