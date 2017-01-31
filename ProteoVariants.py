import sys
sys.path.extend(['/home/walzer/immuno-tools/Fred2'])
__author__ = 'walzer'

from Fred2.IO.ADBAdapter import EIdentifierTypes
import Fred2.Core.Generator as g
from Fred2.IO import FileReader
from Fred2.Core.Allele import Allele
from Fred2.Core.Peptide import Peptide
from Fred2.Core.Peptide import Protein
from Fred2.EpitopePrediction import EpitopePredictorFactory
from Fred2.IO.MartsAdapter import MartsAdapter
from Fred2.Core.Result import EpitopePredictionResult
from Fred2.Core import Variant
from Fred2.Core import VariationType
from Fred2.Core import MutationSyntax
import pandas
import argparse
import os
import logging
import datetime
import itertools
from collections import defaultdict


VERSION = 0.2
PRED_METH = "netmhc"
IDS = dict()
vartype_namemap = {0: "SNP", 1: "DEL", 2: "INS", 3: "FSDEL", 4: "FSINS", 5: "UNKNOWN"}


def peptides_to_fasta(pep_df):
    s = ""
    count = 1
    if type(pep_df) == list:
        for pep in pep_df:
            s += pep_to_fasta_entry(pep, count)
            count += 1
    elif type(pep_df) == EpitopePredictionResult:
        for i, row in pep_df.iterrows():
            pred = ""
            for k in i[1:]:
                pred += k + '=' + str(i[0].get_metadata(k, True)).strip("{}") + ';'
            s += pep_to_fasta_entry(i[0], count, pred)
            count += 1
    else:
        return None
    return s


def pep_to_fasta_entry(pep, count, pred=""):
    annot = ""
    ids = set()
    genes = set()
    for trID, poslist in pep.proteinPos.iteritems():
        ids.add(trID)
        for pos in poslist:
            #using last applied variant only ([-1])
            t = pep.get_variants_by_protein_position(trID, pos).values()[-1][-1].type
            genes.add(pep.get_variants_by_protein_position(trID, pos).values()[-1][-1].gene)
            annot += vartype_namemap[t] + '='
            if t == VariationType.SNP:
                ap = pep.get_variants_by_protein_position(trID, pos).keys()[-1]  # as above - only use [-1] ... this is TODO better
                aseq = str(pep)[:ap] + '|' + str(pep)[ap] + '|' + str(pep)[(ap + 1):] \
                    if ap < len(str(pep)) - 1 \
                    else str(pep)[:ap] + '|' + str(pep)[ap] + '|'
                annot += aseq + "({ms})".format(ms = pep.get_variants_by_protein_position(trID, pos).values()[0][0].coding.values()[0].aaMutationSyntax) + ';'
            else:
                vpos = pep.get_variants_by_protein_position(trID, pos).values()[0][0].coding[trID.split(":")[0]].protPos
                annot += str(vpos) + '><' + str(pos) + ';'
    return "gnl|" + str(count) + "|" + ';'.join(genes) + "|" + ';'.join(ids) + "|" + \
        annot + str(pred) + '\n' + \
        str(pep) + '\n'


def proteins_to_fasta(pro_ls):
    s = ""
    count = 1
    if type(pro_ls) == list:
        for pro in pro_ls:
            s += pro_to_fasta_entry(pro, count)
            count += 1
    elif type(pro_ls) == dict:
        for pId, pro in pro_ls.iteritems():
            s += pro_to_fasta_entry(pro, count)
            count += 1
    else:
        print "What kind of proteins are these?! Abort.", type(pro_ls)
    return s


def pro_to_fasta_entry(pro, count):
    tID = pro.transcript_id.split(":")[0]  # assumes it like NM_014826:FRED2_1 - hell breaks loose if ori. id contains :
    mutsyn = [j.coding[tID].aaMutationSyntax for i in pro.vars.values() for j in i]  # only given transcript - variants contained hold all transcripts mutations syntaxes - Fred2 protein-variant design qwirk
    id = pro.gene_id
    anno = str(pro.get_metadata(PRED_METH, True)).strip("{}")

    #>db|UniqueIdentifier|EntryName ProteinName OS=OrganismName[ GN=GeneName]PE=ProteinExistence SV=SequenceVersion

    return "gnl|" + str(count) + "|" + id + "|" + pro.transcript_id + "=" + '&'.join(mutsyn) + "|" + anno + \
           '\n' + str(pro) + '\n'


def annotate_protein_from_peptides(pep_df):
    if not type(pep_df) == EpitopePredictionResult:
        #print type(pep_df)
        if not ((type(pep_df) == pandas.core.frame.DataFrame) and all([type(i[0]) == Peptide for i, r in pep_df.iterrows()])):
            return None

    pros = dict()
    global PRED_METH
    # alleles
    a = list(set(itertools.chain.from_iterable([r.index for i, r in pep_df.iterrows()])))
    # {allele:{protein_id:[(prot_pos,score), ...]}}
    anno = {allele:defaultdict(list) for allele in a}
    for i, r in pep_df.iterrows():
        pros.update(i[0].proteins)
        PRED_METH = i[1]
        for k, v in i[0].proteinPos.iteritems():
            for idx in r.index:
                anno[idx][k] += [(x, r[idx]) for x in v]


    # invert outer and inner keys of nested dict
    anno_inv = {}
    for k1, subdict in anno.items():
            for k2, v in subdict.items():
                anno_inv.setdefault(k2, {})[k1] = v

    for k, v in anno_inv.iteritems():
        pros[k].log_metadata(PRED_METH, v)

    # all proteins with predictions, metadata = 'prediction method'={allele:[(prot_pos,score), ...]
    return pros


def is_stop_gain(pro):
    tID = pro.transcript_id.split(":")[0]  # assumes it like NM_014826:FRED2_1 - hell breaks loose if ori. id contains :
    mutsyn = [j.coding[tID].aaMutationSyntax for i in pro.vars.values() for j in i ]  # only given transcript - variants contained hold all transcripts mutations syntaxes - Fred2 protein-variant design qwirk
    if any([x.endswith("*") for x in mutsyn]):
        return True
    return False


def __main__():
    parser = argparse.ArgumentParser(version=VERSION)
    parser.add_argument('-V', '--variations', dest="var_file", help='<Required> full path to the input variations', required=True)
    parser.add_argument('-o', "--outfile", dest="outfile_path", help="Created fasta file", required=True)
    parser.add_argument('-d', "--digest", dest="digest", type=int, help="Length of peptides for predigestion and prediction, default 9.")
    parser.add_argument('-a', "--alleles", dest="alleles", help="Input alleles for prediction")
    parser.add_argument('-p', "--predict", dest="predict_with", help="Method of prediction, needs alleles & length, allowed:[{m}]".format(m=PRED_METH))
    parser.add_argument('-f', "--filter", dest="filter", type=float, help="Only include sequences with predictions above the given threshold (e.g. 0.4256 for at least weak binder), needs predict")
    parser.add_argument('-P', "--Proteins", dest="only_proteins", action='store_true', help="Will write only proteins.")
    parser.add_argument('-b', "--base", dest="basefasta_path", help="If given, entries are replaced by the variation.")

    options = parser.parse_args()
    if len(sys.argv) <= 1:
        parser.print_help()
        sys.exit(1)

    if options.filter and not options.predict_with:
        parser.print_help()
        print "Need alleles with predict option, aborting!"
        sys.exit(1)

    if options.predict_with and not options.alleles:
        parser.print_help()
        print "Need alleles with predict option, aborting!"
        sys.exit(1)

    temp_dir = "/tmp/"

    logging.basicConfig(filename=os.path.splitext(options.outfile_path)[0] + "_{:%d-%m-%Y_%H-%M-%S}".format(datetime.datetime.now()) + '.log',
                        filemode='w+', level=logging.DEBUG)  #, format='%(levelname)s:%(message)s'
    logging.info("Starting variant fasta creation " + options.outfile_path + " at " + str(datetime.datetime.now()))
    logging.warning("verbosity turned on")

    #... look at theos filter, ligandoqc, fasta-distributions, lica and the morgenstellen server conten scripts
    # complete proteins?
    # only containing binders?
    # k-mers?
    # binders only?
    # FastaSlicer.py?
    # remove original if homozygous (needs fasta input)?
    # add germline variant option? or expect all to be in one vcf?

# MyObject = type('MyObject', (object,), {})
# options = MyObject()
# setattr(options,"var_file","/home/walzer/immuno-tools/Fred2/Fred2/Data/examples/vcftestfile3.vcf")
#
# vt = os.path.splitext(options.var_file)[-1]
# if ".vcf" == vt:
#     vcfvars, accessions = FileReader.read_vcf(options.var_file)
#
# mart_db = MartsAdapter(biomart="http://grch37.ensembl.org")
#
# transcript_gen = g.generate_transcripts_from_variants(vcfvars, mart_db, id_type=EIdentifierTypes.REFSEQ)
# transcripts = [x for x in transcript_gen if x.vars]
# transcript_gen = g.generate_transcripts_from_variants(vcfvars, mart_db, id_type=EIdentifierTypes.REFSEQ)
# protein_gen = g.generate_proteins_from_transcripts(transcript_gen)
# proteins = [x for x in protein_gen if x.vars]
# for p in proteins:
#     p.gene_id = p.vars.values()[0][0].gene
#
#
# for t in transcripts:
#     t.gene_id = t.vars.values()[0].gene
#

    vt = os.path.splitext(options.var_file)[-1]
    if ".vcf" == vt:
        vcfvars, accessions = FileReader.read_vcf(options.var_file)
    elif ".GSvar" == vt:
        pass
        # vcfvars = FileReader.read_GSvar(options.var_file)
    else:
        m = "Could not read variants {f}, aborting.".format(f=options.var_file)
        logging.error(m)
        print m
        sys.exit(1)

    mart_db = MartsAdapter(biomart="http://grch37.ensembl.org")  # TODO guess id_type for mart_db from accessions

    transcript_gen = g.generate_transcripts_from_variants(vcfvars, mart_db, id_type=EIdentifierTypes.REFSEQ)

    protein_gen = g.generate_proteins_from_transcripts(transcript_gen)
    proteins = [x for x in protein_gen if x.vars]  # removing unvaried

    for p in proteins:
        p.gene_id = p.vars.values()[0][0].gene  # assume gene name from first variant

    proteins = [p for p in proteins if not is_stop_gain(p)]  # kick out stop gains

    # First exit option
    if not (options.predict_with or options.filter) and options.only_proteins:
        if options.basefasta_path:
            # TODO - replace from base fasta
            print "N/A"
            sys.exit(0)
        else:
            e = proteins_to_fasta(proteins)
            with open(options.outfile_path, 'w') as f:
                f.write(e)
            sys.exit(0)

    # From now on, digestion must be set somehow
    if not options.digest:
        digest = 9
    else:
        digest = options.digest
    peptide_gen = g.generate_peptides_from_proteins(proteins, digest)
    peptides = [x for x in peptide_gen]
    peptides_var = [x for x in peptides if any(x.get_variants_by_protein(y) for y in x.proteins.keys())]  # removing unvaried

    # Second exit option
    if not (options.predict_with or options.filter):
        e = peptides_to_fasta(peptides_var)
        with open(options.outfile_path, 'w') as f:
            f.write(e)
        sys.exit(0)

    # From now on, predictions are needed
    try:
        target_alleles_set = set(FileReader.read_lines(options.alleles, in_type=Allele))
    except Exception as e:
        m = "Could not read alleles file {f}, aborting.".format(f=options.alleles)
        logging.error(m)
        print m, "what:", str(e)
        sys.exit(1)

    try:
        ttn = EpitopePredictorFactory(options.predict_with)
    except Exception as e:
        m = "Could not initialize prediction method {f}, aborting.".format(f=options.predict_with)
        logging.error(m)
        print m
        sys.exit(1)

    try:
        preds = ttn.predict(peptides_var, alleles=target_alleles_set)
    except Exception as e:
        print "something went wrong with the prediction", options.inf, options.predict_with, "what:", str(e)
        sys.exit(1)

    # punch prediction results in peptide metadata (inside pandas dataframe)
    #PRED_METH = set()
    for i, row in preds.iterrows():
        for j in i[1:]:
            i[0].log_metadata(j, dict(zip(row.index, row.values)))
            #PRED_METH.add(j)  # need that later

    # Third exit option
    if not options.filter:
        if options.only_proteins:
            if options.basefasta_path:
                # TODO - replace from base fasta plus prediction annotation
                print "N/A"
                sys.exit(0)
            else:
                prs = annotate_protein_from_peptides(preds)
                e = proteins_to_fasta(prs)
                with open(options.outfile_path, 'w') as f:
                    f.write(e)
                sys.exit(0)
        else:
            e = peptides_to_fasta(preds)
            with open(options.outfile_path, 'w') as f:
                f.write(e)
            sys.exit(0)

    # kick out nonbinder
    preds_f = preds[(preds > options.filter).any(axis=1)]

    # Fourth exit option
    if options.only_proteins:
        if options.basefasta_path:
            # TODO - replace from base fasta binders only plus prediction annotation
            print "N/A"
            sys.exit(0)
        else:
            prs = annotate_protein_from_peptides(preds_f)
            e = proteins_to_fasta(prs)
            with open(options.outfile_path, 'w') as f:
                f.write(e)
            sys.exit(0)
    else:
        e = peptides_to_fasta(preds_f)
        with open(options.outfile_path, 'w') as f:
            f.write(e)
        sys.exit(0)

# allele1 = Allele("HLA-A*02:01")
# allele2 = Allele("HLA-A*03:01")
# peptide_gen = g.generate_peptides_from_proteins(proteins, 9)
# peptides = [x for x in peptide_gen]
# peptides_var = [x for x in peptides if any(x.get_variants_by_protein(y) for y in x.proteins.keys())]
# predictor = EpitopePredictorFactory("Syfpeithi")
# results = predictor.predict(peptides_var, alleles=[allele1,allele2])
# results[(results > 25).any(axis=1)]
#
# pep = results[(results > 30).any(axis=1)].iloc[0].name[0]
# #or use
# pep = results[(results > 30).any(axis=1)].index.get_level_values('Seq')[4]
# pep.proteins.keys()
# #['NM_014675:FRED2_0']
#
# results[(results > 30).any(axis=1)].iloc[4].index[0]
# #HLA-A*02:01
# results[(results > 30).any(axis=1)].iloc[4].iloc[0]
# #31
#
# results[(results > 10).any(axis=1)].iloc[0].name[0].log_metadata(results[(results > 10).any(axis=1)].index.get_level_values('Method')[-1], dict(zip(results[(results > 10).any(axis=1)].iloc[0].index, results[(results > 10).any(axis=1)].iloc[0].values)))
#
# for i, row in results[(results > 7).any(axis=1)].iterrows():
#     for j in i[1:]:
#         i[0].log_metadata(j, dict(zip(row.index, row.values)))
#
#
# results[(results > 7).any(axis=1)].iloc[1].name[0].get_metadata("syfpeithi")

# j = None
# for i,row in ppn.iterrows():
#      for trID, poslist in i[0].proteinPos.iteritems():
#              for pos in poslist:
#                      if i[0].get_variants_by_protein_position(trID, pos).values()[-1][-1].type == VariationType.SNP:
#                             print "snp"
#                             j = i

    # def get_decorated_sequence(self):
    #     deco_snp = set()
    #     deco_fs = set()
    #     deco_id = set()
    #     for fid, o in self.proteinPos.iteritems(): # defaultdict(list,{'ENST00000372470:FRED2_0': [508],'ENST00000413998:FRED2_0': [508]})
    #         for p in o:
    #             for vp, vo in self.get_variants_by_protein_position(fid, o[0]).iteritems(): # {514: [Variant(g.43815008G>T)]}
    #                 if vo[0].type == VariationType.SNP: # for now only first variant at this position - todo sort by fs<indel<snp
    #                     deco_snp.add(vp - p)
    #     sorted(deco_snp, reverse=True)
    #     seq = str(self)
    #     for vp in deco_snp:
    #         if vp >= len(str(self)):
    #             logging.warn("unresolvable variant situation")
    #         seq = seq[:vp] + '|' + seq[vp] + '|' + seq[vp + 1:] if vp < len(str(self)) - 1 else seq[:vp] + '|' + seq[vp] + '|'
    #     return seq

# j = None
# for i,row in ppn.iterrows():
#      for trID, poslist in i[0].proteinPos.iteritems():
#              for pos in poslist:
#                      if i[0].get_variants_by_protein_position(trID, pos).values()[-1][-1].type == VariationType.SNP:
#                             print i[0].proteinPos, "snp", i[0].get_variants_by_protein_position(trID, pos).keys(), i[0].get_variants_by_protein_position(trID, pos).values()[0][0].coding.values()[0].aaMutationSyntax
#                             j = i
#
# NM_014675:FRED2_0 2005 snp 3 p.Ser2009Cys
# NM_014675:FRED2_0 2008 snp 0 p.Ser2009Cys
# NM_014675:FRED2_0 2007 snp 1 p.Ser2009Cys
# NM_152888:FRED2_1 689 snp 6 p.Gly696Trp
# NM_152888:FRED2_1 692 snp 3 p.Gly696Trp
# NM_152888:FRED2_1 690 snp 5 p.Gly696Trp
# NM_014675:FRED2_0 2000 snp 8 p.Ser2009Cys
# NM_152888:FRED2_1 693 snp 2 p.Gly696Trp
# NM_014675:FRED2_0 2006 snp 2 p.Ser2009Cys
# NM_152888:FRED2_1 688 snp 7 p.Gly696Trp
# NM_152888:FRED2_1 687 snp 8 p.Gly696Trp
# NM_152888:FRED2_1 691 snp 4 p.Gly696Trp
# NM_014675:FRED2_0 2001 snp 7 p.Ser2009Cys
# NM_014675:FRED2_0 2002 snp 6 p.Ser2009Cys
# NM_152888:FRED2_1 694 snp 1 p.Gly696Trp
# NM_014675:FRED2_0 2004 snp 4 p.Ser2009Cys
# NM_014675:FRED2_0 2003 snp 5 p.Ser2009Cys
# NM_152888:FRED2_1 695 snp 0 p.Gly696Trp


if __name__ == '__main__':
    __main__()

