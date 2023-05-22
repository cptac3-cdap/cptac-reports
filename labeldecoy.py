#!/usr/bin/env python3

import sys
import re
from psmdb.PSMDb import PSMDb

db = PSMDb(filename=sys.argv[1])
decoyprefix = sys.argv[2]
trans=db.begin()
for pep in db.peptides(connection=trans):
    pep.decoy = False
for pr in db.proteins(connection=trans):
    if pr.accession.startswith(decoyprefix):
        pr.decoy=True
        for pep in pr.peptides():
            pep.decoy=True
    else:
        pr.decoy=False
for g in db.genes(connection=trans):
    if g.name.startswith(decoyprefix):
        g.setdata('decoy',True)
    else:
        g.setdata('decoy',False)
for pr in db.proteins(connection=trans):
    if pr.decoy == False:
        for pep in pr.peptides():
            pep.decoy=False
trans.commit(close=True)
