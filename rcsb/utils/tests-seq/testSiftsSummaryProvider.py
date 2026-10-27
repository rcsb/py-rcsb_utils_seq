##
# File:    SiftsSummaryProviderTests.py
# Author:  J. Westbrook
# Date:    11-Dec-2018
#
# Update:
#  20-Oct-2026 dwp Add test for loading a cache built by older code (uppercase entry ID keys)
#  26-Oct-2026 dwp Add test for SIFTS source entry IDs (fails once SIFTS uses extended PDB IDs)
#
#
##
"""
Various utilities for processing SIFTS  correspondence data.

"""

__docformat__ = "google en"
__author__ = "John Westbrook"
__email__ = "jwest@rcsb.rutgers.edu"
__license__ = "Apache 2.0"

import logging
import os
import platform
import shutil
import time
import unittest

from rcsb.utils.io.FileUtil import FileUtil
from rcsb.utils.io.MarshalUtil import MarshalUtil
from rcsb.utils.seq.SiftsSummaryProvider import SiftsSummaryProvider

HERE = os.path.abspath(os.path.dirname(__file__))
TOPDIR = os.path.dirname(os.path.dirname(HERE))

logging.basicConfig(level=logging.INFO, format="%(asctime)s [%(levelname)s]-%(module)s.%(funcName)s: %(message)s")
logger = logging.getLogger()


def getNonShortEntryIds(entryIds):
    """Return the entry IDs that aren't short (4-character) PDB IDs, e.g. extended (12-character) PDB IDs."""
    return sorted({eId for eId in entryIds if len(eId) != 4})


class SiftsSummaryProviderTests(unittest.TestCase):
    def setUp(self):
        self.__dirPath = os.path.join(os.path.dirname(TOPDIR), "rcsb", "mock-data")
        self.__srcDirPath = os.path.join(self.__dirPath, "sifts-summary")
        #
        # self.__cacheDirPath = os.path.join(HERE, "test-output", "sifts-summary")
        self.__cachePath = os.path.join(HERE, "test-output", "CACHE")
        #
        self.__startTime = time.time()
        logger.info("Starting %s at %s", self.id(), time.strftime("%Y %m %d %H:%M:%S", time.localtime()))

    def tearDown(self):
        endTime = time.time()
        logger.info("Completed %s at %s (%.4f seconds)\n", self.id(), time.strftime("%Y %m %d %H:%M:%S", time.localtime()), endTime - self.__startTime)

    def testWriteReadSiftsSummaryCache(self):
        abbreviated = "TEST"
        su = SiftsSummaryProvider(srcDirPath=self.__srcDirPath, cachePath=self.__cachePath, useCache=False, abbreviated=abbreviated)
        eCountW = su.getEntryCount()
        logger.info("SIFTS entry count %d", eCountW)
        self.assertGreaterEqual(eCountW, 140000)
        aL = su.getIdentifiers("1CBS", "A", "UNPID")
        self.assertEqual(len(aL), 1)
        logger.info("Completed %s at %s (%.4f seconds)", self.id(), time.strftime("%Y %m %d %H:%M:%S", time.localtime()), time.time() - self.__startTime)
        su = SiftsSummaryProvider(cachePath=self.__cachePath, useCache=True)
        eCountR = su.getEntryCount()
        logger.info("SIFTS entry count %d", eCountR)
        self.assertEqual(eCountW, eCountR)
        aL = su.getIdentifiers("102M", "A", "UNPAL")
        self.assertEqual(len(aL), 1)
        unpIdL = su.getIdentifiers("102M", "A", "UNPID")
        self.assertEqual(len(unpIdL), 1)
        pfamIdL = su.getIdentifiers("102M", "A", "PFAMID")
        self.assertEqual(len(pfamIdL), 1)
        #
        uniqueUniprotAL = su.getUniqueIdentifiers(idType="UNPID")
        logger.info("Unique UniProt identifiers %d", len(uniqueUniprotAL))
        #
        uniqueUniprotBL = su.getEntryUniqueIdentifiers(su.getEntries(), idType="UNPID")
        self.assertEqual(len(uniqueUniprotAL), len(uniqueUniprotBL))
        #
        saoLD = su.getLongestAlignments("2B7X", ["A", "B", "C", "D"])
        logger.info("2B7X_1 alignments: %r", saoLD)
        #
        if abbreviated == "PROD":
            iproIdL = su.getIdentifiers("102M", "A", "IPROID")
            self.assertEqual(len(iproIdL), 4)
            goIdL = su.getIdentifiers("102M", "A", "GOID")
            self.assertGreaterEqual(len(goIdL), 5)
            logger.debug("GO IDS (%d) %r", len(goIdL), goIdL)
            logger.debug("unpIdl %r pfamIdL %r iprodL %r", unpIdL, pfamIdL, iproIdL)

        logger.info("Completed %s at %s (%.4f seconds)", self.id(), time.strftime("%Y %m %d %H:%M:%S", time.localtime()), time.time() - self.__startTime)

    def testReadLegacySiftsSummaryCache(self):
        """Test that a cache built by older code (uppercase entry ID keys, e.g. restored from the stash) is still usable."""
        cachePath = os.path.join(HERE, "test-output", "CACHE-legacy-sifts")
        shutil.rmtree(cachePath, ignore_errors=True)
        legacyD = {
            "1CBS": {"A": {"UNPID": ["P29373"], "UNPAL": [{"UP": "P29373", "BG": 1, "LEN": 137, "UBG": 1, "UND": 137}], "TAXID": ["9606"]}},
            "1234": {"A": {"UNPID": ["P00001"]}},
        }
        mU = MarshalUtil()
        mU.mkdir(os.path.join(cachePath, "sifts-summary"))
        ok = mU.doExport(os.path.join(cachePath, "sifts-summary", "sifts-summary-py3.pic"), legacyD, fmt="pickle")
        self.assertTrue(ok)
        su = SiftsSummaryProvider(cachePath=cachePath, useCache=True)
        self.assertEqual(su.getEntries(), ["1234", "1cbs"])
        for entryId in ["1CBS", "1cbs"]:
            self.assertEqual(su.getIdentifiers(entryId, "A", "UNPID"), ["P29373"])
            self.assertEqual(su.getTaxIds(entryId, "A"), ["9606"])
            self.assertEqual(su.getAlignmentCount(entryId, "A"), 1)
        self.assertEqual(su.getEntryUniqueIdentifiers(["1234", "1CBS"]), ["P00001", "P29373"])

    @unittest.skipIf(platform.system() != "Darwin", "Skip long development troubleshooting test")
    def testWriteSiftsSummaryCacheJson(self):
        entrySaveLimit = 50
        su = SiftsSummaryProvider(
            srcDirPath=self.__srcDirPath,
            cachePath=self.__cachePath,
            cacheKwargs={"fmt": "json", "indent": 3},
            useCache=False,
            entrySaveLimit=entrySaveLimit,
            abbreviated=False,
        )
        eCount = su.getEntryCount()
        logger.info("SIFTS entry count %d", eCount)
        self.assertGreaterEqual(eCount, entrySaveLimit)

    def testSiftsSourceEntryIds(self):
        """Check that the SIFTS summary source data (current EBI SIFTS files) still uses short (4-character) PDB IDs.

        Production reads the same files from an rsync mirror. Builds the abbreviated ("TEST") summary, as production does,
        from the two files it reads.
        """
        srcDirPath = os.path.join(HERE, "test-output", "sifts-source-data")
        cachePath = os.path.join(HERE, "test-output", "CACHE-sifts-source")
        fU = FileUtil()
        for fn in ["pdb_chain_uniprot.csv.gz", "pdb_chain_pfam.csv.gz"]:
            self.assertTrue(fU.get("https://ftp.ebi.ac.uk/pub/databases/msd/sifts/flatfiles/csv/" + fn, os.path.join(srcDirPath, fn)))
        su = SiftsSummaryProvider(srcDirPath=srcDirPath, cachePath=cachePath, useCache=False, abbreviated="TEST")
        self.assertTrue(su.testCache())
        # Fails once the source data uses extended PDB IDs; then update this package and rcsb.utils.dictionary for them
        badL = getNonShortEntryIds(su.getEntries())
        self.assertEqual(badL, [], "SIFTS source data has entry IDs that aren't 4-char PDB IDs: %r" % badL[:10])


def readSiftsInfo():
    suiteSelect = unittest.TestSuite()
    # suiteSelect.addTest(SiftsSummaryProviderTests("testWriteSiftsSummaryCacheJson"))
    suiteSelect.addTest(SiftsSummaryProviderTests("testWriteReadSiftsSummaryCache"))
    suiteSelect.addTest(SiftsSummaryProviderTests("testReadLegacySiftsSummaryCache"))
    suiteSelect.addTest(SiftsSummaryProviderTests("testSiftsSourceEntryIds"))
    return suiteSelect


if __name__ == "__main__":
    mySuite = readSiftsInfo()
    unittest.TextTestRunner(verbosity=2).run(mySuite)
