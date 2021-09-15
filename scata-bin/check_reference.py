#!/usr/bin/env python

import constants
import sys, time, sge, re, os, traceback


import MySQLdb

from Bio import SeqIO
from Bio.Seq import Seq

from constants import *


refsetid = int(sys.argv[1])

log_entry("Starting %s" % (sys.argv[0]))

log_entry("Connecting to database.")


db = MySQLdb.connect( host=db_host, user=db_user, passwd=db_pass, db=db_db)
db_c = db.cursor(MySQLdb.cursors.DictCursor)


db_c.execute("SELECT owner,Name,ready,CutSequence,CutLength FROM ReferenceSets WHERE idReferenceSet = %s", (refsetid,))

row = db_c.fetchone()
owner = row["owner"]
name = row["Name"]

motif = row["CutSequence"]
bp_to_keep = int(row["CutLength"])


if row["ready"] == 1:
    log_entry("Reference Set %d already set to ready" % (refsetid))
    sys.exit(0)

db_c.execute("SELECT Email FROM Users WHERE idUsers = %s", (owner,))
email = db_c.fetchone()["Email"]


lengths = [ ]
seqs = [ ]
errors = ""
warnings = ""
n=0
m=0

ref_count = 0
mean_len = 0
max_len = 0

try:
    refs = SeqIO.parse(open(refset_dir + "/" + sys.argv[1] + ".fas"),"fasta")
    for seq_record in refs:
        n += 1

        spl = []
        seq = ""
        
        if bp_to_keep > 0:
            spl = re.split(motif.upper(),str(seq_record.seq).upper())

            if len(spl) == 1:
                spl = re.split(motif.upper(),str(seq_record.seq.reverse_complement()).upper())

            if len(spl) == 1:
                continue


            seq = spl[-1][0:bp_to_keep]
            #print seq_record.id, seq
        else:
            seq = str(seq_record.seq).upper()


        m+=1
        l=len(seq)
        lengths.append(l)

        seq_record.seq = Seq(seq, seq_record.seq.alphabet)

        seqs.append(seq_record)
    
        if l > ref_max_len:
            errors += "Sequence %d (%s) is too long (%d > %d).\n" \
                      % (n, seq_record.id, l, ref_max_len)
        elif l > ref_long:
            warnings += "Sequence '%s' (%d, %d bp) is longer than %d bp. Consider cutting down.\n" \
                        % (seq_record.id, n, l, ref_long)
        if len(errors) > 10000:
            errors += "Too many errors, bailing out\n"
            break

        if len(warnings) > 20000:
            errors += "Too many warnings. Please revise your reference set.\n"
            break

        if n > 1000000:
            errors += "More than 1M references. Is this really a relevant reference set?"
            break



        ref_count = len(lengths)
        mean_len = sum(lengths) / ref_count
        max_len = max(lengths)
    os.remove(refset_dir + "/" + sys.argv[1] + ".fas")
    SeqIO.write(seqs, open(refset_dir + "/" + sys.argv[1] + ".fas", "wct"),"fasta")
except:
     exceptionType, exceptionValue, exceptionTraceback = sys.exc_info()
     exception_text = "".join(traceback.format_exception(exceptionType, exceptionValue,
                                    exceptionTraceback))
     open(log_dir + "/check_ref_exceptions_%04d.txt" % (refsetid), "w").write(repr(sys.exc_info()))
     send_mail("mikael.durling@slu.se", "Check RefSet Job %d failed" % (refsetid), 
               "Owner %s\nRefsetid: %d\nException: %s\nBacktrace %s" % (email,
                                                                     refsetid,
                                                                     repr(sys.exc_info()),
                                                                     exception_text))
     errors += "Parse error, is it really a FASTA file?\n"


subject = "SCATA: reference set " + name
body = "The set of references you submitted to SCATA has been processed.\n"
if errors != "":
    log_entry("Errors found in reference set %d for user %s" % (refsetid, email))
    subject += " contains errors"
    body += "Errors were found while processing the reference set. An attempt to give details of the errors is found below. \n\n"
    body += "The reference set has been removed from the database. Pleas correct any errors and retry your upload.\n\n"
    body += "Error log:\n" + errors
    body += "\n\nWarnings:\n" + warnings
    db_c.execute("DELETE FROM ReferenceSets WHERE idReferenceSet = %s", (refsetid,))
    db.commit()
    os.remove(refset_dir + "/" + sys.argv[1] + ".fas")
else:
    log_entry("Reference set %d accepted for user %s. [%d sequences, mean length %d, max length %d]" \
              % (refsetid, email, ref_count, mean_len, max_len))
    subject += " is accepted"
    body += "No errors were found and the reference set can now be used in subsequent analyses.\n\n"
    body += "For your detail, there appears to be %d references with a mean length of %d bp and the longest reference is %d bp.\n" % \
            (ref_count, mean_len, max_len)

    if(bp_to_keep > 0):
        body += "\nIntially there were %d sequences submitted. %d sequences were discarded since they\n" \
                "lack the filtering motif supplied.\n" % (n, (n - m))
    if len(warnings):
        body += "\nHowever, some warnings about the reference set should be noted:\n" + warnings + "\n\n"
    db_c.execute("UPDATE ReferenceSets SET ready=1 WHERE idReferenceSet = %s", (refsetid,))
    db.commit()

send_mail(email,subject,body)
