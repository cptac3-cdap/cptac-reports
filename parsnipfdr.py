#!/bin/env sh

import sys, os, os.path, csv, subprocess, math, re
from operator import itemgetter

paramsfile = sys.argv[1]

psmfile = sys.argv[2]
assert psmfile.endswith('.psm')

params=dict(targetprotfdr=1.0,
            specfdrprecision=2,
            bygene='False',
            parsnip_args="-U 2 -C 1 --tieresolution=Keep --pepweights=FDR",
           )
for l in open(paramsfile):
    key,value = map(str.strip,l.split(None,1))
    params[key] = value

tooldir = os.path.split(os.path.abspath(sys.argv[0]))[0]
# reportsdir = os.path.abspath(os.path.join(tooldir,'..','lib','CPTAC-CDAP-Reports',))
parsnip1 = os.path.abspath(os.path.join(tooldir,'parsnip1'))
mayu =  os.path.abspath(os.path.join(tooldir,'mayu'))
if not os.path.exists(parsnip1):
    parsnip1 = [ "python3", os.path.abspath(os.path.join(tooldir,'parsnip1.py')) ]
    mayu = [ "python3", os.path.abspath(os.path.join(tooldir,'mayu.py')) ]
else:
    parsnip1 = [ parsnip1 ]
    mayu = [ mayu ]

required_params = filter(None,"""
nprot
targetprotfdr
specfdrprecision
bygene
parsnip_args
""".split())

for p in required_params:
    assert p in params

targetprotfdr = float(params['targetprotfdr'])
specfdr = float(params.get('initspecfdr',targetprotfdr))
specfdr_precision = int(params['specfdrprecision'])
specfdr_mult = 10**(specfdr_precision)
specfdr = int(round(specfdr_mult*specfdr))
bygene = eval(params['bygene'])
split1 = re.split(r'(?<![\\])"',params['parsnip_args'])
parsnip_args = []
for i in range(len(split1)):
    if i % 2 == 0:
        parsnip_args.extend(split1[i].split())
    else:
        parsnip_args.append(split1[i])

basename = psmfile.rsplit('.',1)[0]
fdrestfile = basename + ".mayu.tsv"

fdrestheaders = list(map(str.strip,filter(None,"""
  SpecFDR%
  TargetProteins
  DecoyProteins
  RawProtFDR%
  MayuProtFDR%
  Selected
""".split())))

wh = open(fdrestfile,'w')
wh.write("\t".join(fdrestheaders)+"\n")
wh.close()

parsfile = basename + ".pars.txt"

mayuheaders = list(map(str.strip,filter(None,"""
  Filename
  Unshared
  TotalProteins
  TargetProteins
  DecoyProteins
  RawProtFDR%
  MayuProtFDR%
  EstTrueProteins
""".split())))

while True:
    
    # Run parsnip1 then mayu to estimate prot FDR
    args = list(parsnip1)
    args.extend(["-d", psmfile])
    args.extend(["-T", "%.*f"%(specfdr_precision,specfdr/float(specfdr_mult))])
    args.extend(parsnip_args)
    if bygene:
        args.extend(["--bygene"])
    args.extend(["--mode","Stats,Dump"])
    args.extend(["--itersolve","50"])
    args.extend(["--ttotal","15"])
    args.extend(["--noalignments"])
    args.extend(["--extension","pars"])

    print("Execute:"," ".join(args),file=sys.stdout)
    proc = subprocess.Popen(args, stdin=None, shell=False)
    retcode = proc.wait()
    assert(retcode == 0)

    args = list(mayu)
    if bygene:
        args.extend(["--bygene"])
    args.extend(["--accessions",parsfile,"--summary",params['nprot'],psmfile])
    print("Execute:"," ".join(args),file=sys.stdout)
    proc = subprocess.Popen(args, stdin=None, shell=False, stdout=subprocess.PIPE)
    for l in proc.stdout:
        data = dict(zip(mayuheaders,l.decode().split()))
    retcode = proc.wait()
    assert(retcode == 0)
    data['SpecFDR%'] = "%.*f"%(specfdr_precision,specfdr/float(specfdr_mult))
    data['RawProtFDR%'] = data['RawProtFDR%'].rstrip('%')
    data['MayuProtFDR%'] = data['MayuProtFDR%'].rstrip('%')
    data['Selected'] = ""

    wh = open(fdrestfile,'a')
    wh.write("\t".join(map(str,map(data.get,fdrestheaders)))+"\n")
    wh.close()

    values = []
    for r in csv.DictReader(open(fdrestfile),dialect='excel-tab'):
        rvals = list(map(float,map(r.get,("SpecFDR%","MayuProtFDR%"))))
        rvals[0] = int(round(rvals[0]*specfdr_mult))
        values.append(rvals)
    values.sort()

    # Are we done?
    nvalues = len(values)
    specfdrs = list(map(itemgetter(0),values))
    protfdrs = list(map(itemgetter(1),values))

    if max(protfdrs) <= targetprotfdr:
        # We always start with the largest spec FDR we will tolerate...
        chosen_specfdr = max(values,key=lambda t: t[1])[0]
        break

    # figure out the next spec FDR to try...
    if min(protfdrs) > targetprotfdr:
        specfdr = int(round(min(specfdrs)/10.0))
        continue

    maxgoodi = max(filter(lambda i:                  protfdrs[i] <= targetprotfdr, range(nvalues)))

    minbadi =  min(filter(lambda i: i > maxgoodi and protfdrs[i] >  targetprotfdr, range(nvalues)))

    if specfdrs[minbadi]-specfdrs[maxgoodi] <= 1:
        # We are done!
        chosen_specfdr = specfdrs[maxgoodi]
        break

    lb = (specfdrs[maxgoodi],protfdrs[maxgoodi])
    ub = (specfdrs[minbadi], protfdrs[minbadi])
    propdist = min(0.75,max(0.25,(targetprotfdr-lb[1])/(ub[1]-lb[1])))
    newspecfdr = int(round(lb[0]+propdist*(ub[0]-lb[0])))
    if newspecfdr == lb[0]:
        newspecfdr += 1
    elif newspecfdr == ub[0]:
        newspecfdr -= 1
    assert newspecfdr not in specfdrs

    specfdr = newspecfdr

# Re-write the output file with the chosen spec FDR marked...
rows = []
for r in csv.DictReader(open(fdrestfile),dialect='excel-tab'):
    if int(round(specfdr_mult*float(r['SpecFDR%']))) == chosen_specfdr:
        r['Selected'] = 'SELECTED'
    rows.append(r)

wh = open(fdrestfile,'w')
wh.write("\t".join(fdrestheaders)+"\n")
for r in rows:
    wh.write("\t".join(map(r.get,fdrestheaders))+"\n")
wh.close()
