#!/usr/bin/env python
__author__ = 'walzer'
import sys
import os
from datetime import datetime
import logging
import argparse
from os.path import isfile, isdir, join, basename
from os import listdir
import pickle
import pyopenms as oms
from os.path import join
from Fred2.Core.Allele import Allele
from Fred2.Core.Peptide import Peptide
from Fred2.EpitopePrediction import EpitopePredictorFactory

VERSION = "0.2"

#~ path = "/share/projects/ligandomics/BD/MZID/id_sp_xd+PI/"
#~ idf = "150523_DYTK_BD-ZH03_Skin_W_C18postSCX_20%_Rep#1_50cm195min3s_msms1.idXML"
# fin = "/share/projects/ligandomics/BD/MZID/id_sp_xd+PI/150523_DYTK_BD-ZH03_Skin_W_C18postSCX_20%_Rep#1_50cm195min3s_msms1.idXML"

def categorize(score):
    if score >= 0.6384377847127609:
        return "strong"
    elif score >= 0.4256251898085073:
        return "weak"
    else:
        return None

def __main__():
    parser = argparse.ArgumentParser(version=VERSION)
    parser.add_argument('-in',  dest="inf", help='<Required> full path to the input idXML directory', required=True)
    parser.add_argument('-out', dest="out", help="<Required> full path to the output idXML directory", required=True)

    options = parser.parse_args()
    if len(sys.argv) <= 1:
        parser.print_help()
        sys.exit(1)

    if not (options.inf or options.out):
        parser.print_help()
        sys.exit(1)

    if not (isdir(options.inf) or isdir(options.out)):
        parser.print_help()
        sys.exit(1)

    idfs = [join(options.inf, f) for f in listdir(options.inf) if isfile(join(options.inf, f)) and f.endswith('.idXML') and "Tue39L243" not in f]

    pepstr = set()
    
    for inf in idfs:
        pros = list()
        peps = list()
        f = oms.IdXMLFile()
        #~ f.load(join(path, idf), pros, peps)
        f.load(inf, pros, peps)
        for pep in peps:
            for h in pep.getHits():
                if "decoy" not in h.getMetaValue("target_decoy"):
                    if 7 < len(h.getSequence().toUnmodifiedString()) < 12:
                        pepstr.add(h.getSequence().toUnmodifiedString())

    donordir = {"BD-ZH03": ['A*01:01','A*11:01','B*15:01','B*35:01','C*03:03','C*04:01'], 
    "BD-ZH02": ['A*11:01','A*68:01','B*15:01','B*35:03','C*03:03','C*04:01']}

    ttn = EpitopePredictorFactory('netmhc')
    target_alleles_net = None
    for do in donordir:
        if do in options.inf:
            target_alleles_net = donordir[do]

    if not target_alleles_net:
        print "no donor recognized", options.inf
        sys.exit(1)

    hla_pref = 'HLA-'
    target_alleles_net_a = [Allele(hla_pref + x) for x in target_alleles_net]

    es = [Peptide(x) for x in pepstr]

    try:
        preds_n = ttn.predict(es, alleles=target_alleles_net_a)
    except Exception as e:
        print "something went wrong with the netMHC prediction", options.inf, "what:",  e.
        pickle.dump(pepstr, join(options.outf, "seqset.pkl"))
        sys.exit(1)

    preds = dict()    
    for index, row in preds_n.iterrows():
        score = row.max()
        allele = str(row.idxmax())
        categ = categorize(score)
        seq = row.name[0].tostring()
        if categ:
            preds[seq] = (allele,categ,score)

    for inf in idfs:
        pros = list()
        peps = list()
        f = oms.IdXMLFile()
        f.load(inf, pros, peps)
        for pep in peps:
            hits = pep.getHits()
            nhits = list()
            for h in hits:
                if h.getSequence().toUnmodifiedString() in preds:
                    x = preds[h.getSequence().toUnmodifiedString()]
                    h.setMetaValue('binder',x[0])
                    h.setMetaValue(str(x[1]),x[2])
                    nhits.append(h)
                else:
                    nhits.append(h)
            pep.setHits(nhits)                    
        f.store(options.out, pros, peps)


if __name__ == '__main__':
    __main__()