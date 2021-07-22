##
# File:    SeqAlignTests.py
# Author:  j. westbrook
# Date:    19-Sep-2019
# Version: 0.001
#
# Update:
##
"""
Test cases for sequence alignment utility class.

"""

import logging
import os
import unittest

from rcsb.utils.seq.SeqAlign import SeqAlign, splitSeqAlignObjList

HERE = os.path.abspath(os.path.dirname(__file__))
TOPDIR = os.path.dirname(os.path.dirname(os.path.dirname(HERE)))

logging.basicConfig(level=logging.INFO, format="%(asctime)s [%(levelname)s]-%(module)s.%(funcName)s: %(message)s")
logger = logging.getLogger()


class SeqAlignTests(unittest.TestCase):
    def setUp(self):
        self.__exdbL = [
            {"authAsymId": "A", "entitySeqIdBeg": 1, "entitySeqIdEnd": 5, "dbSeqIdBeg": 20, "dbSeqIdEnd": 25, "dbName": "UNP", "dbAccession": "P000001", "entityAlignLength": 5},
            {"authAsymId": "A", "entitySeqIdBeg": 1, "entitySeqIdEnd": 5, "dbSeqIdBeg": 20, "dbSeqIdEnd": 25, "dbName": "UNP", "dbAccession": "P000001", "entityAlignLength": 5},
            {"authAsymId": "A", "entitySeqIdBeg": 1, "entitySeqIdEnd": 100, "dbSeqIdBeg": 20, "dbSeqIdEnd": 120, "dbName": "UNP", "dbAccession": "P000001", "entityAlignLength": 5},
            {
                "authAsymId": "A",
                "entitySeqIdBeg": 200,
                "entitySeqIdEnd": 300,
                "dbSeqIdBeg": 1220,
                "dbSeqIdEnd": 1320,
                "dbName": "UNP",
                "dbAccession": "P000001",
                "entityAlignLength": 5,
            },
            {
                "authAsymId": "A",
                "entitySeqIdBeg": 400,
                "entitySeqIdEnd": 500,
                "dbSeqIdBeg": 1420,
                "dbSeqIdEnd": 1520,
                "dbName": "UNP",
                "dbAccession": "P000001",
                "entityAlignLength": 5,
            },
        ]
        self.__siftsL = [
            {"UP": "P000001", "BG": 1, "UBG": 20, "LEN": 5},
            {"UP": "P000001", "BG": 1, "UBG": 20, "LEN": 5},
            {"UP": "P000001", "BG": 1, "UBG": 100, "LEN": 100},
            {"UP": "P000001", "BG": 50, "UBG": 70, "LEN": 100},
            {"UP": "P000001", "BG": 400, "UBG": 20, "LEN": 50},
            {"UP": "P000001", "BG": 1000, "UBG": 1200, "LEN": 500},
        ]

    def testSeqAlignExDb(self):
        """Test grouping functions for EXDB alignments,"""
        try:
            seqAlignObjL = []
            for alObj in self.__exdbL:
                seqAlignObjL.append(SeqAlign("PDB", **alObj))
            grpD = splitSeqAlignObjList(seqAlignObjL)
            self.assertEqual(len(grpD), 3)
            logger.debug("grpD %r", grpD)
            #
            seqAlignObj = self.__exdbL[0]
            sa = SeqAlign("PDB", **seqAlignObj)
            self.assertEqual(sa.getEntityAlignLength(), 5)
            self.assertEqual(sa.getEntitySeqIdBeg(), 1)
            self.assertEqual(sa.getEntitySeqIdEnd(), 5)
            self.assertEqual(sa.getDbName(), "UNP")
            self.assertEqual(sa.getDbAccession(), "P000001")
            self.assertEqual(sa.getDbSeqIdBeg(), 20)
            self.assertEqual(sa.getDbSeqIdEnd(), 25)
        except Exception as e:
            logger.exception("Failing with %s", str(e))
            self.fail()

    def testSeqAlignSifts(self):
        """Test grouping functions for SIFTS alignments,"""
        try:
            seqAlignObjL = []
            for alObj in self.__siftsL:
                seqAlignObjL.append(SeqAlign("SIFTS", **alObj))
            grpD = splitSeqAlignObjList(seqAlignObjL)
            self.assertEqual(len(grpD), 3)
            logger.info("grpD %r", grpD)
        except Exception as e:
            logger.exception("Failing with %s", str(e))
            self.fail()


def suiteSeqAlignTests():
    suiteSelect = unittest.TestSuite()
    suiteSelect.addTest(SeqAlignTests("testSeqAlignExDb"))
    suiteSelect.addTest(SeqAlignTests("testSeqAlignSifts"))
    return suiteSelect


if __name__ == "__main__":
    mySuite = suiteSeqAlignTests()
    unittest.TextTestRunner(verbosity=2).run(mySuite)
