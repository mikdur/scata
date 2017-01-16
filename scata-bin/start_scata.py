#!/usr/bin/python

import sys, time, sge, os, re, traceback, random
import scata

import MySQLdb

from Bio import SeqIO

from constants import *


class MySQLPrinter:

    def __init__(self,jobid,msg=""):
        self.start_time = time.time()
        sys.stdout.write("Feedback through MySQLPrinter()\n")
        sys.stdout.flush()
        self.last=0
        

    def update_state(self, now, total):
        sys.stdout.write("\r" + (" " * self.last) + "\r")
        runtime = time.time() - self.start_time
        tpu = runtime / (float(now)+0.000000000001)
        to_go = tpu * (total - now)
        eta = ""
        try:
            eta = time.ctime(time.time() + to_go)
            if to_go > (30 * 24 * 3600):
                eta = "Estimating ETA"
        except ValueError:
            eta = "Unknown"
        msg = self.msg + ("ETA: " + eta)
        self.last = len(msg)
        sys.stdout.write(msg)
        sys.stdout.flush()
        db = MySQLdb.connect( host=db_host, user=db_user, passwd=db_pass, db=db_db)
        db_c = db.cursor(MySQLdb.cursors.DictCursor)
	db_c.execute("SELECT * from Jobs WHERE idJobs = %s", (jobid,))
	row = db_c.fetchone()
	if not row:
	    sys.stdout.write("Job gone, exiting")
	    exit(0)
        db_c.execute("UPDATE Jobs SET description=%s WHERE idJobs = %s", (msg,jobid,))
        db.commit()

    def set_msg(self, msg):
        self.msg=msg

    def reset(self):
        self.start_time = time.time()

jobid = int(sys.argv[1])

log_entry("Connecting to database.")

log_text = ""
db = MySQLdb.connect( host=db_host, user=db_user, passwd=db_pass, db=db_db)
db_c = db.cursor(MySQLdb.cursors.DictCursor)


db_c.execute("SELECT owner,name,status,locked,tagSet,jobSettingsSet,sgeProject  FROM Jobs JOIN Users on (Jobs.owner = Users.idUsers) WHERE idJobs = %s", (jobid,))

row = db_c.fetchone()
owner = row["owner"]
name = row["name"]
sge_proj = row["sgeProject"]

if not (row["status"] == 1 and row["locked"] == 1):
    log_entry("Dataset %d not in state to run" % (jobid))
    sys.exit(0)

db_c.execute("SELECT Email FROM Users WHERE idUsers = %s", (owner,))
email = db_c.fetchone()["Email"]

log_entry("Preparing job %d for user %s (%d)" % (jobid, email, owner))

# Parameters that the is not allowed to touch:

config = { "job_id" :         { "value": "scata%04d" % (jobid),
                                "descr": "Job ID to identify job" },
           "job_name" :       { "value": name.replace(": ", ":"),
                                "descr": "Job name" },
           "sge_params" :     { "value": sge_params_scata_jobs + sge_proj,
                                "descr": "Per job parameters to SGE" },
           "graph":           { "value": "0",
                                "descr": "Create GV files" },
           "work_dir" :       { "value": run_dir + "/run_" + str(jobid),
                                "descr": "Working directory where temporary working files are stored" },
           "lib_dir":         { "value": base_dir,
                                "descr": "Path to directory with all scripts" },
           "output_dir":      { "value": "%s/scata%04d" % (run_dir, jobid),
                                "descr": "Directory where result files are stored" },
           "tags":            { "value": tagset_dir + "/" + str(row["tagSet"]) + ".txt",
                                "descr": "File with tagset definition" },
           "formatdb":        { "value": formatdb_path,
                                "descr": "Path to formatdb executable" },
           "usearch":         { "value": usearch_path,
                                "descr": "Path to usearch" },
           "vsearch":         { "value": vsearch_path,
                                "descr": "Path to usearch" },
           "blastall":        { "value": blastall_path,
                                "descr": "Path to blastall executable" },
           "muscle":          { "value": "muscle -diags -maxiters 2",
                                "descr": "Path and options to muscle"},
           "polyh_num":       { "value": "3",
                                "descr": "Deprecated option of polyh scoring" },
           "polyh_weight":    { "value": "1",
                                "descr": "Deprecated option of polyh scoring" }
           }



# Add datasets to config
db_c.execute("""
SELECT Name, idDatasets FROM Datasets d
   JOIN JobDatasets jd
      ON (d.idDatasets = jd.dataSetId)
   WHERE jd.jobId = %s""", (jobid,) );

names = [ ]
d_ids = [ ]
for dataset in db_c:
    names.append(re.sub("[^-_.,A-Za-z0-9]", "_", dataset["Name"]))
    d_ids.append(str(dataset["idDatasets"]))

d_seq = map(lambda a: "%s/%s.seqs.pick" % (dataset_dir, a), d_ids)
d_stat = map(lambda a: "%s/%s.stat.pick" % (dataset_dir, a), d_ids)


config["454names"] = { "value":  " ".join(names),
                       "descr":  "Names of datasets" }
config["454seq"] = { "value": " ".join(d_seq),
                       "descr": "Seqs pick files with sequence data" }
config["454stat"] = { "value": " ".join(d_stat),
                      "descr": "Seq stats pick files with stats" }

# Add reference sets to config
db_c.execute("SELECT refId FROM JobRefs WHERE jobId = %s", (jobid,) );

r_ids = [ ]
for refset in db_c:
    r_ids.append(str(refset["refId"]))

r_fas = map(lambda a: "%s/%s.fas" % (refset_dir, a), r_ids)


config["reference_seqs"] = { "value":  " ".join(r_fas)  if len(r_fas) else "none",
                             "descr":  "Paths to reference sequences" }
                            

# Pull out all user options from the tables.
db_c.execute("""
SELECT p.Name, IF(p.type='text',psv.textValue,pv.value) AS value, p.keyword, p.coerce, p.type
	FROM ParameterSets ps 
	JOIN ParameterSettoValue psv 
		ON (ps.idSettingsSet = psv.parameterSetId) 
	JOIN ParameterValue pv 
		ON (psv.parameterValueId = pv.idValue) 
	JOIN Parameters p 
		ON (pv.parameterId = p.idParameter) 
	WHERE ps.idSettingsSet=%s""", (row["jobSettingsSet"],) )

for setting in db_c:
    if setting["type"] == "text":
        try:
            if setting["coerce"] == "float":
                foo = float(setting["value"])
            elif setting["coerce"] == "int":
                foo = int(setting["value"])
        except ValueError:
            log_text += "Setting '%s' is not of correct format. ('%s' can't be treatad as %s)" % \
                        (setting["Name"], setting["value"], setting["coerce"])
            log_entry("Setting '%s' is not of correct format. ('%s' can't be treatad as %s)" % 
                        (setting["Name"], setting["value"], setting["coerce"]))
    config[setting["keyword"]] = {"value": setting["value"],
                                 "descr": setting["Name"] }


log_entry("Creating config file")
config_file = open(run_dir + "/config_%04d.txt" % (jobid),"w")

for option, value in config.iteritems():
    config_file.write("# %s\n%s: %s\n\n" % (value["descr"], option, value["value"]))
config_file.close()

# Call SCATA proper to run the clustering process
if len(log_text) == 0:
    try:
        db_c.execute("UPDATE Jobs SET status=2, description='Preparing job',startTime=NOW() WHERE idJobs = %s", (jobid,))
        db.commit()
        config_file = run_dir + "/config_%04d.txt" % (jobid)


        # Set up printer
        pr = MySQLPrinter(jobid)
        # Run scata
        scata.run_scata(config_file,pr)

        # Reconnect to db as a lot of time might have passed
        db = MySQLdb.connect( host=db_host, user=db_user, passwd=db_pass, db=db_db)
        db_c = db.cursor(MySQLdb.cursors.DictCursor)

        # zip result file and copy to results dir.
        filename = "%s_" % (config["job_id"]["value"]) + \
                   "".join(random.sample('abcdefghijklmnopqrstuvxyzABCDEFGHIJKLMNOPQRSTUVXYZ',10)) + \
                   ".zip"

        full_filename = result_dir + "/" + filename
        res_dir = "scata%04d" % (jobid)
        os.system("cd %s && zip -r %s %s" % ( run_dir, full_filename, res_dir))
        log_entry("Job %d done" % (jobid))
        send_mail(email, "Your job %s is ready" % (name), "Your job has been successfully processed. " \
                  "The results are available for download on the SCATA website.\n")
        db_c.execute("UPDATE Jobs SET status=3, filename=%s, endTime=NOW(), description=%s WHERE idJobs = %s",
                     (filename, ("Done (" + time.ctime() + ")"), jobid))
        os.system("rm -rf %s %s %s" % (config["work_dir"]["value"], config["output_dir"]["value"], config_file) )
        db.commit()
    except SystemExit:
        send_mail("mikael.durling@slu.se", "Job %d deleted while running" % (jobid), 
                  "Owner %s\nJobId: %d\n" % (email,
                                             jobid))
        os.system("rm -rf %s %s %s" % (config["work_dir"]["value"], config["output_dir"]["value"], config_file) )
    except:
        exceptionType, exceptionValue, exceptionTraceback = sys.exc_info()
        exception_text = "".join(traceback.format_exception(exceptionType, exceptionValue,
                                          exceptionTraceback))
        open(log_dir + "/exceptions_%04d.txt" % (jobid), "w").write("%s\n\n%s" % (repr(sys.exc_info()),
                                                                                  exception_text))
        send_mail("mikael.durling@slu.se", "Job %d failed" % (jobid), 
                  "Owner %s\nJobId: %d\nException: %s\nBacktrace %s" % (email,
                                                         jobid,
                                                         repr(sys.exc_info()),
                                                         exception_text))
        log_entry("Job %d for %s failed" % (jobid, email))
        log_text+="An error occurred during your run. The administrator has been notified.\n"
        log_text+="Please try again later, this can be a temporary error. Contact the administrator \n"
        log_text+="if there are any questions."
        send_mail(email, "Job %s failed" % (name), log_text)
        db_c.execute("UPDATE Jobs SET status=4 WHERE idJobs = %s", (jobid,))
        db.commit()

else:
    send_mail(email, "Job %d could not be started." % (jobid),
              "Your job could not be started. There were errors in the parameter set you " + \
              "selected. Please create a new parameter set with the errors corrected.\n\nLog:\n" + \
              log_text)
    log_entry("Job %d could not be started." % (jobid))
    db_c.execute("UPDATE Jobs SET status=4 WHERE idJobs = %s", (jobid,))
    db.commit()
    exit(0)



log_entry("Done.")



        
        
