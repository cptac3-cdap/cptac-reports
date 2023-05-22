#!/usr/bin/env python3

import sys, os, csv, os.path, re
from collections import defaultdict

def sortkey(tag):
    m = re.search(r'^(\w+)-(\d+)(N|C)?$',tag)
    assert m, repr(tag)
    return (int(m.group(2)),-1*(m.group(3)=="")+1*(m.group(3)=="C"))

def intstrkey(val):
    try:
        return (int(val),"")
    except (ValueError,TypeError) as e:
        pass
    return (1e+20,val)

# labels file/database
labels = defaultdict(dict)
if len(sys.argv) >= 3 and os.path.exists(sys.argv[2]):
    labellot = None
    for i,row in enumerate(csv.reader(open(sys.argv[2]),dialect='excel-tab')):
        if len(row) < 5:
            continue
        if row[0].startswith(">"):
            labellot = row[0][1:].strip()
            headers = list(row)
            headers[0] = "label"
            continue
        elif row[0].strip() == "" and i == 0:
            labellot = os.path.split(sys.argv[2])[1].rsplit('.',1)[0]
            headers = list(row)
            headers[0] = "label"
            continue
        data = dict(list(zip(headers,row)))
        labels[labellot][data['label']] = data

# Sample/experimental design file
rows = list(csv.reader(open(sys.argv[1]),dialect='excel-tab'))
headers = rows[0]
rows = rows[1:]

if 'AnalyticalSampleIndex' not in headers:
    asampord = 1
    asamporddict = {}
    for i,r in enumerate(rows):
        d = dict(list(zip(headers,r)))
        if d.get('AnalyticalSample') != None:
            asamp = d.get('AnalyticalSample')
            if asamp not in asamporddict:
                asamporddict[asamp] = asampord
                asampord += 1
            r.append(asamporddict[asamp])
    headers.append('AnalyticalSampleIndex')

# Files to analyze
anyratios = False
out = []
freq = defaultdict(dict)
prefix = None
alltags = []
for l in sorted(sys.stdin):
    filename = os.path.split(l)[1]
    basename,mzid,gz = filename.rsplit('.',2)
    filematched = False
    for i,r in enumerate(rows):
        d = dict(list(zip(headers,r)))
        if re.search(d['FileNameRegEx'],basename):
            filematched = True
            labellot = d.get('LabelReagent')
            ratios = defaultdict(list)
            if 'Ratios' in d:
                for j,(numerator,denominator) in enumerate([tuple(r.split('/')) for r in d.get('Ratios',"").split(',')]):
                    numerator = numerator.strip()
                    denominator = denominator.strip()
                    if denominator not in ratios[numerator]:
                        ratios[numerator].append((j+1,denominator))
            if labellot and labellot not in labels:
                labels[labellot] = dict()
                for i1,row1 in enumerate(csv.reader(open(labellot+'.txt'),dialect='excel-tab')):
                    if i1 == 0:
                        headers1 = list(row1); headers1[0] = 'label'
                        continue
                    d1 = dict(list(zip(headers1,row1)))
                    labels[labellot][d1['label']] = d1
            if prefix == None:
                tags = set([k for k in headers if re.search(r"^\d+",k)])
                if '114' in tags:
                    if len(tags) == 4:
                        prefix = "iTRAQ"
                    elif len(tags) == 8:
                        prefix = "iTRAQ8"
                elif '126' in tags or '126C' in tags:
                    if len(tags) == 6:
                        prefix = "TMT6"
                    elif len(tags) == 10:
                        prefix = "TMT10"
                    elif len(tags) == 11:
                        prefix = "TMT11"
                    elif len(tags) == 16:
                        prefix = "TMT16"
                    elif len(tags) == 18:
                        prefix = "TMT18"
                    else:
                        raise RuntimeError("Bad tag!?!?!")
                elif '000' in tags:
                    prefix = "NoTag"
                else:
                    raise RuntimeError("Bad tag!?!?!")
                for tag in tags:
                    tag = "%s-%s"%(prefix,tag)
                    alltags.append(tag)
                alltags.sort(key=sortkey)
            if d.get('AnalyticalSample'):
                asamp = d.get('AnalyticalSample')
                asampord = int(d.get('AnalyticalSampleIndex'))
            else:
                asamp = tuple([d.get(tag.split('-')[1]) for tag in alltags])
                asampord = 0
            if d.get('Fraction'):
                fraction = d.get('Fraction')
            else:
                fracregex = d.get('FractionRegex',r'\_[fF]?([A-Z]|[0-9]+)$')
                m = re.search(fracregex,basename)
                if not m:
                    raise RuntimeError("Can't extract fraction number from basename: %s"%(basename,))
                fraction = m.group(1)
                try:
                    fraction = str(int(fraction))
                except ValueError:
                    pass
            for j,tag in enumerate(alltags):
                k = tag.split('-')[1]
                row = [basename,asamp,asampord,fraction,tag,d[k]]
                if labellot:
                    percs = labels[labellot][k]
                    simplecount = 0
                    simplemap = {"-2": "-2:13C+13C",
                                 "-1": "-1:13C",
                                 "+1": "+1:13C",
                                 "+2": "+2:13C+13C",
                                 "1":  "+1:13C",
                                 "2":  "+2:13C+13C"}
                    for of in simplemap.keys():
                        if of in percs:
                            simplecount += 1
                            percs[simplemap[of]] = percs[of]
                            del percs[of]
                    if simplecount not in (0,4):
                        raise RuntimeError('Header "-2", "-1", "+1", or "+2" missing from labeling reagent file for lot %s'%(offset,labellot,))
                    elif simplecount == 4:
                        percs["-2:13C+15N"] = 0.0
                        percs["-1:15N"] = 0.0
                        percs["+1:15N"] = 0.0
                        percs["+2:15N+13C"] = 0.0
                      
                    for offset in ("-2:13C+13C","-2:13C+15N","-1:13C","-1:15N","+1:15N","+1:13C","+2:15N+13C","+2:13C+13C"):
                        if offset not in percs:
                            raise RuntimeError("Header '%s' missing from labeling reagent file for lot %s"%(offset,labellot,))
                        row.append(float(percs[offset]))

                if k in ratios:
                    denoms = []
                    for ordinal,denom in ratios[k]:
                        try:
                            fdenom = float(denom)
                        except ValueError:
                            fdenom = None
                        if fdenom != 1.0:
                            denoms.append("%d:%s-%s"%(ordinal,prefix,denom))
                        else:
                            denoms.append("%d:1.0"%(ordinal,))
                    row.append(";".join(denoms))
                    anyratios = True
                else:
                    row.extend("")
                out.append(row)
                if (asamp,tag) not in freq[d[k]]:
                    freq[d[k]][(asamp,tag)] = (i,j)
            break
    assert filematched, "File %s not matched to an analytical sample."%(filename,)

sampmap = dict()
for s in freq:
    if len(freq[s]) > 1 and s not in ("POOL","COMMON_INTERNAL_REFERENCE","COMMON_IR","CIR","COMMON_REFERENCE","COMMON_REF"):
        for i,(asamp,tag) in enumerate(sorted(freq[s],key=lambda k: freq[s][k])):
            sampmap[(asamp,tag)] = "%s.%d"%(s,i+1)
    else:
        for asamp,tag in freq[s]:
            sampmap[(asamp,tag)] = s

out.sort(key=lambda t: (t[2],t[1],intstrkey(t[3]),t[0],sortkey(t[4])))
headers = ["spectrafile","analyticalsample","analyticalsampleordinal","fraction"]
if len(alltags) > 1:
    headers.extend(["label","sample"])
else:
    headers.extend(["sample"])
if len(labels) > 0:
    headers.extend(["-2:13C+13C","-2:13C+15N","-1:13C","-1:15N","+1:15N","+1:13C","+2:15N+13C","+2:13C+13C"])
if anyratios:
    headers.extend(["denominator"])
print(",".join(headers))
for r in out:
    r[5] = sampmap[(r[1],r[4])]
    if len(alltags) == 1:
        del r[4]
    if not isinstance(r[1],str):
        r[1] = ":".join([sampmap[(r[1],t)] for t in alltags])
    print(",".join(map(str,r)))
