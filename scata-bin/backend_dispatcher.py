#!/usr/bin/env python

import constants
import sys, time, sge

import MySQLdb



from constants import *



log_entry("Starting %s" % (sys.argv[0]))


log_entry("Entering main loop")

cpu = 0
wallclock = 0

last_job_num = 0


class SgeWrapper:
    def __init__(self):
        self.jobs = [ ]
        self.cpu = 0
        self.wallclock = 0
        self.last_jobnum = 0


    def run(self, cmd, args, name, sge_par, rt, mem):
        sge_job = sge.SGEJob(name, job_script_dir, sge_par, rt, mem)
        sge_job.add_task(cmd, args)
        self.jobs.append(dict(job = sge_job,
                              name = name))
        sge_job.start()
        log_entry("Sent off '%s' to SGE" % (name))

    def step(self):
        if len(self.jobs):
            if len(self.jobs) != self.last_jobnum:
                log_entry("%d jobs are waiting in SGE" % len(self.jobs))
                self.last_jobnum = len(self.jobs)
            
        to_delete = []
        for i in range(len(self.jobs)):
            if self.jobs[i]["job"].check()["done"]:
                self.cpu += self.jobs[i]["job"].get_time()["cpu"]
                self.wallclock += self.jobs[i]["job"].get_time()["wallclock"]
                #self.jobs[i]["job"].clean()
                if self.jobs[i]["job"].get_num_failed() > 0:
                    log_entry("Job '%s' failed, retrying" % (self.jobs[i]["name"]))
                    self.jobs[i]["job"].reset_failed()
                    self.jobs[i]["job"].start()
                else:
                    to_delete.append(i)
                    log_entry("Job '%s' done." % (self.jobs[i]["name"]))
                log_entry("Current resource usage: CPU: %f, Wallclock: %f" % (self.cpu, self.wallclock))
        to_delete.reverse()
        for d in to_delete:
            del self.jobs[d]

sge_wrapper = SgeWrapper()

while True:

    #log_entry("Connecting to database.")


    db = MySQLdb.connect( host=db_host, user=db_user, passwd=db_pass, db=db_db)
    db_c = db.cursor(MySQLdb.cursors.DictCursor)

    n = 0
    db_c.execute("SELECT idTagSets, sgeProject from TagSets JOIN Users on (TagSets.owner = Users.idUsers) WHERE ready=0 AND locked=0")
    for row in db_c:
        sge_wrapper.run(base_dir + "/check_tagset.py", [ str(row["idTagSets"]) ],
                        "ScataT%d" % (row["idTagSets"]), sge_params_backend + str(row["sgeProject"]), 60, "1G")
        n += 1
        db_c.execute("UPDATE TagSets SET locked=1 WHERE idTagSets = %s", (row["idTagSets"],))
        db.commit()

    if n:
        log_entry("%d tagsets to check" % (n))


    db_c.execute("SELECT idReferenceSet, sgeProject  from ReferenceSets JOIN Users on (ReferenceSets.owner = Users.idUsers) WHERE ready=0 AND locked=0")
    n = 0
    for row in db_c:
        sge_wrapper.run(base_dir + "/check_reference.py", [ str(row["idReferenceSet"]) ],
                        "ScataR%d" % (row["idReferenceSet"]), sge_params_backend + str(row["sgeProject"]), 20 * 60, "5G")
        n += 1
        db_c.execute("UPDATE ReferenceSets SET locked=1 WHERE idReferenceSet = %s", (row["idReferenceSet"],))
        db.commit()

    if n:
        log_entry("%d reference to check" % (n))


    #log_entry("Checking data sets")
    db_c.execute("SELECT idDatasets, sgeProject  from Datasets JOIN Users on (Datasets.owner = Users.idUsers) WHERE ready=0 and locked=0")
    n = 0
    for row in db_c:
        sge_wrapper.run(base_dir + "/check_dataset.py", [ str(row["idDatasets"]) ],
                        "ScataD%d" % (row["idDatasets"]), "-R y -V " + sge_params_backend + str(row["sgeProject"]),
#			 "ScataD%d" % (row["idDatasets"]), "-V " + sge_params_backend + str(row["sgeProject"]),
                        96 * 60, "10G")
        n += 1
        db_c.execute("UPDATE Datasets SET locked=1 WHERE idDatasets = %s", (row["idDatasets"],))
        db.commit()

    if n: 
        log_entry("%d datasets to check" % (n))


    #log_entry("Checking for new jobs")
    db_c.execute("SELECT idJobs, sgeProject  from Jobs JOIN Users on (Jobs.owner = Users.idUsers) WHERE status=1 AND locked=0")
    n = 0
    for row in db_c:
        sge_wrapper.run(base_dir + "/start_scata.py", [ str(row["idJobs"])],
                        "ScataJ%d" % (row["idJobs"]), sge_params_scata + str(row["sgeProject"]), 60 * 24 * 20, "10G")
        db_c.execute("UPDATE Jobs SET locked=1 WHERE idJobs = %s", (row["idJobs"],))
        db.commit()
        n += 1
    

    # Step the SGE wrapper one step
    sge_wrapper.step()


    time.sleep(60)
    
