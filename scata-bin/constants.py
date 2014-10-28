# File defining all system constants.


db_host = "scata.mykopat.slu.se"
db_user = "scata_dev"
db_pass = ""
db_db   = "scata_dev"


base_dir = "/proj/mykopat-scata/dev/scata/scata-bin"
tagset_dir = "/mykopat/scata/dev/tagsets"
refset_dir = "/mykopat/scata/dev/referencesets"
dataset_dir = "/mykopat/scata/dev/datasets"
file_dir = "/mykopat/scata/dev/files"
result_dir = "/mykopat/scata/dev/results"
run_dir = "/mykopat/scata/dev/run"
log_dir = "/mykopat/scata/dev/log"

ref_long = 1200
ref_max_len=2000

job_script_dir = "/mykopat/scata/dev/tmp"


sge_params_backend = '-V -R y -w n -p -500 -P '
sge_params_scata_jobs = '-V -R n -w n -r y -p -600 -P '
sge_params_scata = '-V -R n -w n -p -700 -q scata@* -P '

mail_from = "SCATA <noreply@mykopat.slu.se>"

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
    s = smtplib.SMTP("my-mgrid.grid.mykopat.slu.se")
    s.sendmail("www-data@heterobasidion.mykopat.slu.se", [to], msg.as_string())
    s.quit()
    
