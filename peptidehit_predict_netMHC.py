#!/usr/bin/env python
__author__ = 'walzer'
import sys
import os
from datetime import datetime
import logging
import argparse
import pickle
import pyopenms as oms
from os.path import join
from Fred2.Core.Allele import Allele
from Fred2.Core.Peptide import Peptide
from Fred2.EpitopePrediction import EpitopePredictorFactory

VERSION = "0.2"

def categorize(score):
    if score >= 0.6384377847127609:
        return "strong"
    elif score >= 0.4256251898085073:
        return "weak"
    else:
        return None

def __main__():
    parser = argparse.ArgumentParser(version=VERSION)
    parser.add_argument('-in',  dest="inf", help='<Required> full path to the input idXML', required=True)
    parser.add_argument('-out', dest="out", help="<Required> full path to the output idXML", required=True)

    options = parser.parse_args()
    if len(sys.argv) <= 1:
        parser.print_help()
        sys.exit(1)

    if not (options.inf or options.out):
        parser.print_help()
        sys.exit(1)
    #~ else:
    #~ print dir(options)

    #~ if options.temp_path:
        #~ temp_dir = options.temp_path
    #~ else:
        #~ temp_dir = "/tmp/"

    if "Tue39L243" in options.inf:
        print "class2 run recognized", options.inf
        sys.exit(1)

    donordir = {"BD-ZH03": ['A*01:01','A*11:01','B*15:01','B*35:01','C*03:03','C*04:01'], 
    "BD-ZH02": ['A*11:01','A*68:01','B*15:01','B*35:03','C*03:03','C*04:01'],
    "BD-ZH01": ['A*02:01','A*11:01','B*27:05','B*51:01','C*01:02','C*15:02'], 
    "BD-ZH06": ['A*03:01','A*68:02','B*07:02','B*14:02','C*07:02','C*08:02']}

    ttn = EpitopePredictorFactory('netmhc')
    target_alleles_net = None
    for do in donordir:
        if do in options.inf:
            target_alleles_net = donordir[do]

    if not target_alleles_net:
        print "no donor recognized", options.inf
        sys.exit(1)

    pros = list()
    peps = list()
    f = oms.IdXMLFile()
    #~ f.load(join(path, idf), pros, peps)
    f.load(options.inf, pros, peps)

    pepstr = set()
    for pep in peps:
        for h in pep.getHits():
            if "decoy" not in h.getMetaValue("target_decoy"):
                ps = h.getSequence().toUnmodifiedString()
                if 7 < len(ps) < 12:
                    if "U" not in ps:
                        pepstr.add(h.getSequence().toUnmodifiedString())

    hla_pref = 'HLA-'
    target_alleles_net_a = [Allele(hla_pref + x) for x in target_alleles_net]

    es = [Peptide(x) for x in pepstr]

    try:
        preds_n = ttn.predict(es, alleles=target_alleles_net_a)
    except Exception as e:
        print "something went wrong with the netMHC prediction", options.inf, "what:",  str(e)
        pickle.dump(pepstr, open( '.'.join(options.out.split(".")[:-1]) + ".pkl", "wb" ) )
        sys.exit(1)

    preds = dict()    
    for index, row in preds_n.iterrows():
        score = row.max()
        allele = str(row.idxmax())
        categ = categorize(score)
        seq = row.name[0].tostring()
        if categ:
            preds[seq] = (allele,categ,score)

    npeps = list()
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

    #~ f.store(join(path, odf), pros, peps)
    f.store(options.out, pros, peps)


if __name__ == '__main__':
    __main__()