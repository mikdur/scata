# File defining all system constants.


db_host = ""
db_user = "scata"
db_pass = ""
db_db   = "scata2"


base_dir = "/mykopat/scata/scata-bin"
tagset_dir = "/mykopat/scata/tagsets"
refset_dir = "/mykopat/scata/referencesets"
dataset_dir = "/mykopat/scata/datasets"
result_dir = "/mykopat/scata/results"
run_dir = "/mykopat/scata/run"
log_dir = "/mykopat/scata/log"

ref_long = 1200
ref_max_len=2000

job_script_dir = "/mykopat/scata/tmp"


sge_params_backend = '-V -R y -w n -p -500 -P '
sge_params_scata_jobs = '-V -R n -w n -r y -p -600 -P '
sge_params_scata = '-V -R n -w n -p -700 -q scata@* -P '

mail_from = "SCATA <noreply@mykopat.slu.se>"

import sys, time

sys.path.append("/mykopat/Linux-x86_64/lib/python2.5/site-packages/")
sys.path.append("/mykopat/Darwin-x86_64/lib/python2.5/site-packages/")

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
    
