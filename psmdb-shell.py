#!/bin/env python27
import sys, os
from psmdb.PSMDb import PSMDb
from psmdb.model import *
from sqlobject import *
db=PSMDb(filename=sys.argv[1],memory=False,processConnection=True)
if len(sys.argv) >= 3:
    exec(compile(open(sys.argv[2]).read(), sys.argv[2], 'exec'))
    sys.exit(0)
import code
import psmdb.history
code.interact(local=locals())
