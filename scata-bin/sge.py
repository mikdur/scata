import copy, plistlib, subprocess, re, time, sys, os, stat, pickle
from xml.etree import ElementTree
from  xml.parsers.expat import ExpatError


task_template = """#!/bin/sh

#$ -l h_vmem=%s,h_rt=0:%d:0
#$ -wd /mykopat/scata/tmp
#$ -o %s
#$ -N '%s'
#$ -M mikael.durling@mykopat.slu.se
#$ -m a
#$ -V
#$ -j y
#$ %s
%s

. /opt/Modules/default/init/bash
module load mpi4py biopython
cd /mykopat/scata/tmp


%s
"""


class Printer:

    def __init__(self,msg=""):
        self.start_time = time.time()
        self.last = 0
        self.msg = msg
        #sys.stdout.write("Feedback through Printer()\n")
        #sys.stdout.flush()
        

    def update_state(self, now, total):
        sys.stdout.write("\r" + (" " * self.last) + "\r")
        runtime = time.time() - self.start_time
        tpu = runtime / (float(now)+0.000000000001)
        to_go = tpu * (total - now)
        eta = ""
        try:
            eta = time.ctime(time.time() + to_go)
        except ValueError:
            eta = "Unknown"
        msg += ("ETA: " + eta)
        self.last = len(msg)
        sys.stdout.write(msg)
        sys.stdout.flush()
        
    def set_msg(self, msg):
        self.msg=msg

    def reset_timer(self):
        self.start_time = time.time()

        
class SGEJob:
    

    def __init__(self, descr,wd, params="",runtime=2,vf="2G",printer=None):
        self.descr = descr
        self.wd = wd
        self.prefix = wd + "/" + self.descr.replace("/","_").replace(":","_").replace(" ","_").replace("_mykopat_scata_run_config_","") + "_" + str(time.time())
        self.tasks = [ ]
        self.tasks_by_id = {}
        self.locked = False
        self.job_id = 0
        self.task_id = 0
        self.sge_task_ids = { }
        self.sge_id_2_id = {}
        self.check_run = 1
        self.files = []
        self.params = params
        self.array_job = True
        self.no_array=False
        self.runtime = runtime
        self.vf = vf

        if printer == None:
            self.printer = Printer()
        else:
            self.printer = printer
        
        self.sge_path = "/opt/sw/sge/8.1.6/bin/lx-amd64/"

    def add_task(self, cmd, args, dep=[]):
        if self.locked: raise Exception("Can't add tasks to locked job")
        self.task_id += 1
        self.tasks.append({ "cmd" : cmd,
                            "args" : copy.copy(args),
                            "dep" : dep,
                            "id" : str(self.task_id),
                            "status" : "NS",
                            "sge_task_id" : 0,
                            "check" : 0,
                            "done" : False,
                            "wallclock" : 0.0,
                            "cpu" : 0.0,
                            "resets" : 0},
                          )
        self.tasks_by_id[self.tasks[-1]["id"]] = self.tasks[-1]
        return str(self.task_id)

    def start(self):

        if len(self.tasks) == 0:
            return
        if self.locked:
            print "Already started"
            return
        self.locked=True

        self.sge_job = 0

        # Check if there are any dependencies, otherwise make it an array job
        if self.no_array:
            self.array_job=False
        else:
            for t in self.tasks:
                if len(t["dep"]):
                    self.array_job=False
                    break

        if len(self.tasks) == 1:
            self.array_job=False
        
        if self.array_job:
            return self.start_array_job()

        for t in self.tasks:
            if t["done"]:
                continue
            
            deps = ""
            if len(t["dep"]):
                deps = ",".join(map(lambda a: self.sge_task_ids[a], t["dep"]))
                deps = "#$ -hold_jid " + deps
            #print t
            task_script = task_template % ( self.vf,
                                            self.runtime,
                                            self.prefix + "_sge$JOB_ID." + t["id"] + ".out",
                                            self.descr.replace("/","_").replace(":","_").replace("_mykopat_scata_run_config_","") + "_id_" + t["id"],
                                            self.params,
                                            deps,
                                            t["cmd"] + " '" + "' '".join(t["args"]) + "'")
                                            

            self.files.append(self.prefix + "_" + t["id"] + ".out")
            script = open(self.prefix + "_" + t["id"] + "." + "script.sh","w")
            self.files.append(self.prefix +  "_" + t["id"] + "." + "script.sh")
            script.write(task_script)
            script.close()
            os.chmod(self.prefix + "_" + t["id"] + "." + "script.sh", stat.S_IWUSR | stat.S_IRUSR | stat.S_IXUSR)
            restart_count = 0
            while True:
		print self.sge_path + "qsub " + self.prefix + "_" + t["id"] + "." + "script.sh"
                pipe = subprocess.Popen(self.sge_path + "qsub " + self.prefix + "_" + t["id"] + "." + "script.sh" , 
                                        shell=True, bufsize=1024, 
                                        stdout=subprocess.PIPE).stdout
                for line in pipe:
                    #print line
                    if line.find("job") >= 0:
                        t["sge_task_id"] = line.split()[2]
                        self.sge_task_ids[t["id"]]=t["sge_task_id"]
                        self.sge_id_2_id[t["sge_task_id"]]=t["id"]
                t["status"]="S"
                if t["sge_task_id"]:
                    break
                else:
                    if restart_count > 3:
                        raise Exception("Too many restart retries. Cluster down?")
                    print "Retrying to start job"
                    restart_count += 1
                    
        #print self.sge_task_ids, self.tasks

    def start_array_job(self):
        # Create a dictionary to translate between task id and parameters for task
        tasks = {}
        for t in self.tasks:
            tasks[t["id"]] = { "cmd" : t["cmd"] ,
                               "args" : t["args"],
                               "done" : t["done"]}
            
        pickle.dump(tasks, open(self.prefix + "tasks.pick", "wct"))

        # Create runner python script
        runner = """#!/usr/bin/python
import pickle, sys, os
tasks = pickle.load(open('%s'))
t=tasks[sys.argv[1]]

if t['done']:
	exit(0)
        
print 'About to run: ' +t['cmd'] + " '" + "' '".join(t["args"]) + "'"
sys.stdout.flush()
exit((os.system(t['cmd'] + " '" + "' '".join(t["args"]) + "'") & 0xff00) >> 8)


        """ % (self.prefix + "tasks.pick")
        run_script = open(self.prefix + "_runner.py","wct")
        run_script.write(runner)
        run_script.close()
        os.chmod(self.prefix + "_runner.py",  stat.S_IWUSR | stat.S_IRUSR | stat.S_IXUSR)
        self.files.append(self.prefix + "_runner.py")


        # Create task_list
        t_list = "#$ -t 1-" + str(len(self.tasks))
        

        task_script = task_template % ( self.vf,
                                        self.runtime,
                                        self.prefix + "_sge$JOB_ID.$TASK_ID.out",
                                        self.descr.replace("/","_").replace(":","_").replace("_mykopat_scata_run_config_",""),
                                        self.params,
                                        t_list,
                                        self.prefix + "_runner.py $SGE_TASK_ID")
                                            

        self.files.append(self.prefix + ".out")
        script = open(self.prefix +  "_script.sh","w")
        self.files.append(self.prefix +  "_script.sh")
        script.write(task_script)
        script.close()
        os.chmod(self.prefix + "_script.sh", stat.S_IWUSR | stat.S_IRUSR | stat.S_IXUSR)
        restart_count=0
        while True:
            pipe = subprocess.Popen(self.sge_path + "qsub " + self.prefix + "_script.sh" , 
                                    shell=True, bufsize=1024, 
                                    stdout=subprocess.PIPE).stdout
            for line in pipe:
                #print line
                if line.find("job") >= 0:
                    self.sge_job = line.split()[2].split(".")[0]
                    #print self.sge_job
            if self.sge_job:
                break
            else:
                if restart_count > 3:
                    raise Exception("Too many restart retries. Cluster down?")
                print "Retrying to start job"
                restart_count += 1
                    
        
        
    def check(self):
        if len(self.tasks) == 0:
            return {"done": True,
                    "wallclock": 0.0,
                    "cpu": 0.0,
                    "states": { }}
        
        if not self.locked:
            raise Exception("Cannot check() an unlocked job")
        
        if self.array_job:
            return self.check_array_job()

        self.check_run += 1

        pipe = subprocess.Popen(self.sge_path + "qstat -xml -t", 
                                shell=True, bufsize=1024, 
                                stdout=subprocess.PIPE).stdout

        etr = ElementTree.parse(pipe).getroot()

        for e in etr.find("queue_info").findall("job_list"):
            job_id = e.find("JB_job_number").text
            if job_id in self.sge_id_2_id:
                job_id = self.sge_id_2_id[job_id]
                self.tasks_by_id[job_id]["check"] = self.check_run
                self.tasks_by_id[job_id]["status"] = e.find("state").text
            
        
        for e in etr.find("job_info").findall("job_list"):
            job_id = e.find("JB_job_number").text
            if job_id in self.sge_id_2_id:
                job_id = self.sge_id_2_id[job_id]
                self.tasks_by_id[job_id]["check"] = self.check_run
                self.tasks_by_id[job_id]["status"] = e.find("state").text
            
        for t in filter(lambda a: not a["done"] and a["check"] < self.check_run, self.tasks):


            pipe = subprocess.Popen(self.sge_path + "qacct -j " + t["sge_task_id"] + " 2>/dev/null", 
                                shell=True, bufsize=1024, 
                                stdout=subprocess.PIPE).stdout
            for line in pipe:
                l = line.split()
                if l[0] == "ru_wallclock":
                    t["wallclock"] = float(l[1].replace("s",""))
                    t["done"]=True
                    t["status"]="done"
                if l[0] == "cpu":
                    t["cpu"] = float(l[1].replace("s",""))
                    t["done"]=True
                    t["status"]="done"
                if l[0] == "exit_status":
                    t["exit_status"] = int(l[1])
                    t["done"]=True
                    t["status"]="done"
                    
        summary = { "done" : False,
                    "wallclock" : 0.0,
                    "cpu" : 0.0,
                    "states" : { }}
        done = 0
        for t in self.tasks:
            summary["wallclock"] += t["wallclock"]
            summary["cpu"] += t["cpu"]
            if t["status"] in summary["states"]:
                summary["states"][t["status"]] += 1
            else:
                summary["states"][t["status"]] = 1
        
        if "done" in summary["states"] and summary["states"]["done"] == len(self.tasks):
            summary["done"] = True
        return summary

    def check_array_job(self):
        if not self.locked:
            raise Exception("Cannot check() an unlocked job")
        
        self.check_run += 1

        pipe = subprocess.Popen(self.sge_path + "qstat -xml -t", 
                                shell=True, bufsize=1024, 
                                stdout=subprocess.PIPE).stdout

        etr = ElementTree.parse(pipe).getroot()

        for e in etr.find("queue_info").findall("job_list"):
            job_id = e.find("JB_job_number").text
            if job_id == self.sge_job:
                task = e.find("tasks").text
                self.tasks_by_id[task]["check"] = self.check_run
                self.tasks_by_id[task]["status"] = e.find("state").text
            

        for e in etr.find("job_info").findall("job_list"):
            job_id = e.find("JB_job_number").text
            if job_id == self.sge_job:
                tasks = []
                for i in e.find("tasks").text.split(","):
                    if len(i.split(":")) == 2:
                        t = i.split(":")[0].split("-")
                        tasks += [x for x in range(int(t[0]),int(t[1]))]
                    else:
                        tasks.append(int(i))
                    
                state = e.find("state").text
                for id in tasks:
                    self.tasks_by_id[str(id)]["check"] = self.check_run
                    self.tasks_by_id[str(id)]["status"] = e.find("state").text
            
        pipe = subprocess.Popen(self.sge_path + "qacct -j " + self.sge_job + " 2>/dev/null", 
                                shell=True, bufsize=1024, 
                                stdout=subprocess.PIPE).stdout
        for line in pipe:
            l = line.split()
            if l[0] == "taskid":
                t = self.tasks_by_id[l[1]]

            if l[0] == "ru_wallclock":
                t["wallclock"] = float(l[1].replace("s",""))
                t["done"]=True
                t["status"]="done"
            if l[0] == "cpu":
                t["cpu"] = float(l[1].replace("s",""))
                t["done"]=True
                t["status"]="done"
            if l[0] == "exit_status":
                t["exit_status"] = int(l[1])
                t["done"]=True
                t["status"]="done"
            
        summary = { "done" : False,
                    "wallclock" : 0.0,
                    "cpu" : 0.0,
                    "states" : { }}
        done = 0
        for t in self.tasks:
            summary["wallclock"] += t["wallclock"]
            summary["cpu"] += t["cpu"]
            if t["status"] in summary["states"]:
                summary["states"][t["status"]] += 1
            else:
                summary["states"][t["status"]] = 1
        
        if "done" in summary["states"] and summary["states"]["done"] == len(self.tasks):
            summary["done"] = True
        return summary

    
    def wait(self, do_print=True):
        stat = self.check()
        msg = ""
        while not stat["done"]:
            if do_print:
                s = stat["states"]
                num_done = 0
                num_tot = 0
                if "done" in s:
                    num_done += s["done"]
                    num_tot += s["done"]
                if "qw" in s:
                    num_tot += s["qw"]
                if "r" in s:
                    num_tot += s["r"]
                self.printer.update_state(num_done, num_tot)


            while True:
                time.sleep(20)
                try:
                    stat=self.check()
                except ExpatError:
                    if do_print:
                        sys.stdout.write(".")
                        sys.stdout.flush()
                    stat=0
                if stat != 0:
                    break
                

        msg = "\nDone. Wallclock: %.3f cpu: %.3f\r" % ( stat["wallclock"],
                                                        stat["cpu"])
        if do_print:
            sys.stdout.write(msg)
            sys.stdout.write("\n")
            sys.stdout.flush()


    def get_num_failed(self):
        stat = self.check()
        if not stat["done"]:
            return -1

        num_failed = 0
        for t in self.tasks:
            if t["exit_status"] != 0:
                num_failed += 1
        return num_failed

    def reset_failed(self):
        self.locked = False
        self.no_array = True

        increase_runtime = False
        for t in self.tasks:
            t["resets"] += 1
            if t["resets"] >= 6:
                raise Exception("Task " + str(self.sge_job) + " reset too many times.\n" + repr(t))
            
            if t["exit_status"] != 0:
                t["done"]=False
                
                if t["exit_status"] == 137:
                    increase_runtime=True

        if increase_runtime:
            print "some tasks need longer runtime, increasing"
            self.runtime *= 2
                    
        
    def get_time(self):
        time.sleep(20)
        return {"wallclock" : self.check()["wallclock"],
                "cpu" : self.check()["cpu"] }


    def clean(self):
        for f in self.files:
            try:
                os.remove(f)
            except OSError:
                continue
                
