##
# File:    SiftsSummaryProviderTests.py
# Author:  J. Westbrook
# Date:    11-Dec-2018
#
# Update:
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
import time
import unittest

from rcsb.utils.seq.SiftsSummaryProvider import SiftsSummaryProvider

HERE = os.path.abspath(os.path.dirname(__file__))
TOPDIR = os.path.dirname(os.path.dirname(HERE))

logging.basicConfig(level=logging.INFO, format="%(asctime)s [%(levelname)s]-%(module)s.%(funcName)s: %(message)s")
logger = logging.getLogger()


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


def readSiftsInfo():
    suiteSelect = unittest.TestSuite()
    # suiteSelect.addTest(SiftsSummaryProviderTests("testWriteSiftsSummaryCacheJson"))
    suiteSelect.addTest(SiftsSummaryProviderTests("testWriteReadSiftsSummaryCache"))
    return suiteSelect


if __name__ == "__main__":
    mySuite = readSiftsInfo()
    unittest.TextTestRunner(verbosity=2).run(mySuite)
