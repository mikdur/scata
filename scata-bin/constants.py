# File defining all system constants.
import sys

db_host = ""
db_user = ""
db_pass = ""
db_db   = ""


base_dir = "/scata/scata-system/scata-bin"
tagset_dir = "/scata/scata-data/tagsets"
refset_dir = "/scata/scata-data/referencesets"
dataset_dir = "/scata/scata-data/datasets"
result_dir = "/scata/scata-data/results"
file_dir = "/scata/scata-files/files"
run_dir = "/scata/scata-run/run"
log_dir = "/scata/scata-run/log"

# Paths
formatdb_path = "/scata/scata-system/bin/formatdb"
usearch_path = "/scata/scata-system/bin/usearch"
vsearch_path = "/scata/scata-system/bin/vsearch"
blastall_path = "/scata/scata-system/bin/blastall"

ref_long = 3000
ref_max_len=6000

job_script_dir = "/scata/scata-run/tmp"


sge_params_backend = '-V -R y -w n -p -500 -P '
sge_params_scata_jobs = '-V -R n -w n -r y -p -600 -P '
sge_params_scata = '-V -R n -w n -p -700 -q scata@* -P '

mail_from = "SCATA <noreply@slu.se>"

import sys, time


def log_entry(text):
    log_file = open(log_dir + "/scata.log", "a")
    log_file.write("%s: %s %s: %s\n" % (time.ctime(time.time()), sys.argv[0],
                                          sys.argv[1] if len(sys.argv) == 2 else "", text))


import smtplib
from email.mime.text import MIMEText


def send_mail(to, subject, body):
    msg = MIMEText(body,"plain","iso-8859-1")
    msg["To"] = to
    msg["Subject"] = subject
    msg["From"] = mail_from

    #sys.stdout.write(msg.as_string())
    s = smtplib.SMTP("my-gridfront.grid.mykopat.slu.se")
    s.sendmail("www-data@my-scata.mykopat.slu.se", [to], msg.as_string())
    s.quit()
    
