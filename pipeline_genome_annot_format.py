#!/usr/bin/env python
from pipeline_common import *;
C = init_conf()


gffout = C.__genome_annot_format_output__();

cmd = 'gffread ' \
    + '-E %s ' % C.genome_annot \
    + '-o %s'  % gffout;

retval = run_cmd(cmd);

if retval != 0:
    print "ERROR: gffread fails, simply copying genome annotation file instead!!!!!"
    cmd = 'cp %s %s' % (C.genome_annot, gffout)
    retval = run_cmd(cmd)

sys.exit(retval);
