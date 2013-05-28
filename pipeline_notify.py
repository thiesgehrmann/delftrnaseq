#!/usr/bin/python

import os;
import sys;
import smtplib;
from email.mime.text import MIMEText;

from pipeline_common import *;

###############################################################################

def usage(a1):
  print "Usage:  %s <config file> <step> <log file> <status>" % a1;
#edef

if len(os.sys.argv) != 5:
  usage(os.sys.argv[0]);
  os.sys.exit(1);
#fi

###############################################################################

C = PIPELINECONF(os.sys.argv[1]);
step = os.sys.argv[2];
logf = os.sys.argv[3];
stat = os.sys.argv[4];
#stat = "COMPLETED" if stat == "0" else "FAILED";

try:
  msg = MIMEText("Step [%s] in '%s' has %s. You can find the standard output in %s.\n" % (step, C.jobname, stat, logf));
  msg['Subject'] = "%s: [%s] %s" % (C.jobname, step, stat);
  msg['From']    = C.email;
  msg['To']      = C.email;

  s = smtplib.SMTP(C.__pipeline_mail_server__);
  s.starttls();
  s.login(C.__pipeline_mail_user__, C.__pipeline_mail_pass__);
  s.sendmail(C.__pipeline_email__, C.email, msg.as_string());
  s.quit();
  sys.exit(0);
except NameError:
  sys.exit(0);
#etry

