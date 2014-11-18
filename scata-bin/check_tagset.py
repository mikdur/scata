#!/usr/bin/env python
import constants
import sys, time, sge, re, os


import MySQLdb

from constants import *


tagsetid = int(sys.argv[1])

log_entry("Starting %s" % (sys.argv[0]))

log_entry("Connecting to database.")


db = MySQLdb.connect( host=db_host, user=db_user, passwd=db_pass, db=db_db)
db_c = db.cursor(MySQLdb.cursors.DictCursor)


db_c.execute("SELECT owner,tagSetName,ready FROM TagSets WHERE idTagSets = %s", (tagsetid,))

row = db_c.fetchone()
owner = row["owner"]
name = row["tagSetName"]

if row["ready"] == 1:
    log_entry("Tagset %d already set to ready" % (tagsetid))
    sys.exit(0)

db_c.execute("SELECT Email FROM Users WHERE idUsers = %s", (owner,))
email = db_c.fetchone()["Email"]


tag_len = 0
tags = { }

errors = ""
try:
    infile = open(tagset_dir + "/" + sys.argv[1] + ".txt")

    line_no = 0
    for line in infile:
        line_no += 1

        if line == "":
            continue
        if line[-1] == "\n":
            line = line[:-1]
        if line!="" and line[-1] == "\r":
            line = line[:-1]


        ll = line.split(";")

        if len(ll) < 2:
            continue

        if tag_len == 0:
            tag_len = len(ll[1])

        if len(ll[1]) != tag_len:
            errors += "Tag on line %d is not %d bp, as tags on previous lines\n" % (line_no, tag_len)

        ll[1] = ll[1].upper()
        
        if not re.match("^[ACTG]+$",ll[1]):
            errors += "Tag on line %d contains invalid bases\n" % (line_no)

        if not re.match("^[-A-Za-z0-9_]+$",ll[0]):
            errors += "Tag name on line %d contains invalid characters\n" % (line_no)

        if len(ll) > 2:
            for i in range(2, len(ll)):
                if not re.match("^[-A-Za-z0-9_]+$",ll[i])
                    errors += "Pairing tag name %d on line %d contains invalid characters\n" % (i + 1, line_no)
                    
        if len(errors) > 2000:
            errors += "Too many errors, bailing out\n"
            break

        tags[ll[0]]=ll[1]
except:
    errors += "Undefined error, bailing out\n"

subject = "SCATA: tagset " + name
body = "The set of tags you submitted to SCATA has been processed.\n"
if errors != "":
    log_entry("Errors found in tagset %d for user %s" % (tagsetid, email))
    subject += " contains errors"
    body += "Errors were found while processing the tag set. An attempt to give details of the errors is found below.\n\n"
    body += "The tag set has been removed from the database. Pleas correct any errors and retry your upload.\n\n"
    body += "Error log:\n" + errors
    db_c.execute("DELETE FROM TagSets WHERE idTagSets = %s", (tagsetid,))
    db.commit()
    os.remove(tagset_dir + "/" + sys.argv[1] + ".txt")
else:
    log_entry("Tagset %d accepted for user %s. [%d tags, %d bp each]" % (tagsetid, email, len(tags), tag_len))
    subject += " is accepted"
    body += "No errors were found and the tag set can now be used in subsequent analyses.\n\n"
    body += "For your detail, there appears to be %d tags of %d bp each in the tag set.\n" % (len(tags), tag_len)
    db_c.execute("UPDATE TagSets SET ready=1 WHERE idTagSets = %s", (tagsetid,))
    db.commit()

send_mail(email,subject,body)
