
import sqlite3, re
from .SearchResult import *

query = """
SELECT psm.id,
       psm.qvalue as qvalue,
       psm.rank as rank,
       psm.observedneutralmass as precursorMass,
       psm.monoisotopicmasserror as massdiff,
       substr( protdata.Sequence, pep2prot.[OFFSET] + 1, pep2prot.length ) AS pepseq,
       psm.charge as charge,
       spec.precursormz as precursorMz,
       spec.nativeID as nativeID,
       source.name as basename,
       group_concat( DISTINCT ( prot.accession ) ) as proteins,
       group_concat( DISTINCT ( pepmod.site ||( pepmod.[OFFSET] + 1 ) || ':' || mod.MonoMassDelta )  ) as mods
  FROM PeptideSpectrumMatch AS psm,
       Peptide AS pep,
       PeptideInstance AS pep2prot,
       Protein AS prot,
       ProteinData AS protdata,
       Spectrum AS spec,
       SpectrumSource AS source
       LEFT JOIN PeptideModification AS pepmod
              ON pepmod.peptidespectrummatch = psm.id
       LEFT JOIN Modification AS mod
              ON pepmod.modification = mod.id
 WHERE pep.id = psm.peptide
       AND
       pep.id = pep2prot.peptide
       AND
       prot.id = pep2prot.protein
       AND
       psm.spectrum = spec.id
       AND
       spec.source = source.id
       AND
       protdata.id = prot.id
       AND
       prot.isdecoy = 0
       AND psm.id >= ? and psm.id < ?
 GROUP BY psm.id
 ORDER BY spec.source, spec.id
       """

def IDPDB(filename):
    conn = sqlite3.connect(filename)
    conn.row_factory = sqlite3.Row
    conn.text_factory = str
    cur = conn.cursor()
    lastspec = None
    result = None
    blockstart = 0
    blocksize = 1000
    any = True
    while any:
        any = False
        for row in cur.execute(query,(blockstart,blockstart+blocksize)):
            any = True
            if (row['basename'],row['nativeID']) != lastspec:
                if result != None:
                    result.set('ranking_key','qvalue')
                    result.sortPSMs()
                    yield result
                lastspec = (row['basename'],row['nativeID'])
                result = SearchResult()
                pmz = float(row['precursorMz'])
                result.set('precursorMz', pmz)
                result.set('basename', row['basename'])
                m = re.search(r'scan=(\d+)',row['nativeID'])
                assert m
                result.set('start_scan', int(m.group(1)))
                result.set('end_scan', int(m.group(1)))
                result.set('spectrum','%s.%d.%d'%(row['basename'],int(m.group(1)),int(m.group(1))))
            psm = PeptideSpectrumMatch()
            psm.set('peptide',row['pepseq'])
            z = int(row['charge'])
            psm.set('assumed_charge',z)
            psm.set('precursor_neutral_mass', float(row['precursorMass']))
            psm.set('calc_neutral_pep_mass', float(row['precursorMass'])-float(row['massdiff']))
            mods = row['mods']
            if mods == None:
                psm.set('mods','-')
            else:
                psm.set('mods', mods)
            psm.set('protein',','.join([a.rstrip('|') for a in row['proteins'].split(',')]))
            psm.set('qvalue',float(row['qvalue']))
            psm.set('rank',int(row['rank']))
            result.addPSM(psm)
        if result != None:
            result.set('ranking_key','qvalue')
            result.sortPSMs()
            yield result
        blockstart += blocksize

if __name__ == '__main__':

    import sys

    for sr in IDPDB(sys.argv[1]):
        print(sr)
