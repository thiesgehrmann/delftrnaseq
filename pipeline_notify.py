#!/usr/bin/python

import os;
import sys;
import smtplib;
from email.mime.text import MIMEText;

from pipeline_common import *;

###############################################################################

def usage(a1):
  print "Usage:  %s <config file> <step> <log file> <status> <start_time> <end_time>" % a1;
#edef

if len(os.sys.argv) != 7:
  usage(os.sys.argv[0]);
  os.sys.exit(1);
#fi

###############################################################################

C = PIPELINECONF(os.sys.argv[1]);
step = os.sys.argv[2];
logf = os.sys.argv[3];
stat = os.sys.argv[4];
d_s  = os.sys.argv[5];
d_e  = os.sys.argv[6];

try:
  log = open(logf);
  logtext = [ log.readline() for i in xrange(20) ];
  log.close();
  msg = MIMEText("Step [%s] in <%s> has %s.\n\nStarted at: %s\nEnded at: %s\n\nYou can find the standard output in %s.\nFirst 20 lines:\n\n%s\n" % (step, C.jobname, stat, d_s, d_e, logf, '\n'.join(logtext)));
  msg['Subject'] = "DELFT RNA-SEQ PIPELINE: <%s> update" % (C.jobname);
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

