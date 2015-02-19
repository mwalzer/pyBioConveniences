from unittest import TestCase
import warnings
import FastaDB
import tempfile
# import hashlib
import filecmp

__author__ = 'walzer'


class TestFastaDB(TestCase):
    def test_readSeqs(self):
        fdb = FastaDB.FastaDB()

        # case empty
        with warnings.catch_warnings(record=True) as w:
            #  all warnings to be triggered.
            warnings.simplefilter("always")
            #  trigger warning.
            fdb.read_seqs('')
            # tests
            self.assertFalse(len(w) < 1, 'More warnings then expected')
            self.assertTrue(issubclass(w[-1].category, UserWarning))
            self.assertTrue("Could not read file" in str(w[-1].message))
            w.pop()

        # case three.fasta
            fdb = FastaDB.FastaDB()
            fdb.read_seqs('three.fasta')
            self.assertTrue(len(fdb.collection) == 3)
            self.assertTrue(len(fdb.idx) == 3+1)
            self.assertTrue(len(fdb.accs) == 3)
            ts = 0
            for i in fdb.collection.values():
                ts += len(i)
            ts += (len(fdb.collection) - 1)
            self.assertTrue(len(fdb.searchstring) == ts)

        #case addition of three.dat
            fdb.read_seqs('three.dat')
            self.assertTrue(len(fdb.collection) == 6)
            self.assertTrue(len(fdb.idx) == 6+1)
            self.assertTrue(len(fdb.accs) == 6)
            ts = 0
            for i in fdb.collection.values():
                ts += len(i)
            ts += (len(fdb.collection) - 1)
            self.assertTrue(len(fdb.searchstring) == ts)

    def test_writeSeqs(self):
        fdb = FastaDB.FastaDB()
        fdb.read_seqs('three.fasta')
        with tempfile.NamedTemporaryFile() as temp:
            fdb.write_seqs(temp.name)
            # self.assertTrue(hashlib.sha256(open(temp.name, 'rb').read()).digest() == hashlib.sha256(open("three.fasta", 'rb').read()).digest())
            filecmp.cmp(temp.name, 'three.fasta')


    def test_exists(self):
        fdb = FastaDB.FastaDB()
        fdb.read_seqs('three.fasta')
        exi = "TERNEKKQQMGKEYREKIEAEL"
        ine = "XXXNEKKQQMGKEYREKIEAEL"
        self.assertTrue(fdb.exists(exi))
        self.assertFalse(fdb.exists(ine))

    def test_search(self):
        fdb = FastaDB.FastaDB()
        fdb.read_seqs('three.fasta')
        eal = "ELD"
        self.assertDictEqual(fdb.search(eal), {'ELD': 'sp|P31946|1433B_HUMAN'})

    def test_search_all(self):
        fdb = FastaDB.FastaDB()
        fdb.read_seqs('three.fasta')
        eal = "ELD"
        self.assertDictEqual(fdb.search_all(eal), {'ELD': 'sp|P31946|1433B_HUMAN,sp|Q04917|1433F_HUMAN,sp|P62258|1433E_HUMAN'})
