##
# File:    UniProtIdMappingProviderProviderTests.py
# Author:  J. Westbrook
# Date:    30-Nov-2020
#
# Update:
#   10-Dec-2020 jdw add tests for tdd cache persistence
#
##
"""
Tests utilities to manage access to UniProtIdMappingProvider -

These tests from scratch require a couple of hours to run.
(HP 12 core Intel(R) Xeon(R) CPU E5-2680 v3 @ 2.50GHz)

For taxonomy mapping alone -

For cached data performance the "tdd" format data loads fastest (294s + 56 GB mem)
                                "pickle" format data loads (324s + 69 GB mem)

Building cached files (neglecting download of raw mapping files)

   tdd     1937s + 56 GB mem
   pickle  1380s + 69 GB mem

"""

__docformat__ = "google en"
__author__ = "John Westbrook"
__email__ = "jwest@rcsb.rutgers.edu"
__license__ = "Apache 2.0"

import logging
import os
import platform
import resource
import time
import unittest

from rcsb.utils.config.ConfigUtil import ConfigUtil
from rcsb.utils.seq.UniProtIdMappingProvider import UniProtIdMappingProvider

HERE = os.path.abspath(os.path.dirname(__file__))
TOPDIR = os.path.dirname(os.path.dirname(HERE))

logging.basicConfig(level=logging.INFO, format="%(asctime)s [%(levelname)s]-%(module)s.%(funcName)s: %(message)s")
logger = logging.getLogger()


class UniProtIdMappingProviderTests(unittest.TestCase):
    skipFull = True

    def setUp(self):
        self.__cachePath = os.path.join(HERE, "test-output")
        #
        configPath = os.path.join(HERE, "test-data", "stash-config-example.yml")
        self.__configName = "site_info_configuration"
        self.__cfgOb = ConfigUtil(configPath=configPath, defaultSectionName=self.__configName)
        #
        self.__startTime = time.time()
        logger.info("Starting %s at %s", self.id(), time.strftime("%Y %m %d %H:%M:%S", time.localtime()))

    def tearDown(self):
        unitS = "MB" if platform.system() == "Darwin" else "GB"
        rusageMax = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
        logger.info("Maximum resident memory size %.4f %s", rusageMax / 10 ** 6, unitS)
        endTime = time.time()
        logger.info("Completed %s at %s (%.4f seconds)", self.id(), time.strftime("%Y %m %d %H:%M:%S", time.localtime()), endTime - self.__startTime)

    @unittest.skipIf(skipFull, "Very long test")
    def testUniProtIdMappingProviderCachePic(self):
        umP = UniProtIdMappingProvider(self.__cachePath)
        umP.reload(useCache=True, useLegacy=False, fmt="pickel")
        logger.info("Completed %s at %s (%.4f seconds)", self.id(), time.strftime("%Y %m %d %H:%M:%S", time.localtime()), time.time() - self.__startTime)
        #
        taxId = umP.getMappedId("Q6GZX0", mapName="NCBI-taxon")
        logger.info("TaxId %r", taxId)
        self.assertEqual(taxId, "654924")
        #
        taxId = umP.getMappedIdLegacy("Q6GZX0", mapName="NCBI-taxon")
        logger.info("TaxId %r", taxId)
        # self.assertEqual(taxId, "654924")

    @unittest.skipIf(skipFull, "Very long test")
    def testUniProtIdMappingProviderCacheTdd(self):
        umP = UniProtIdMappingProvider(self.__cachePath)
        umP.reload(useCache=True, useLegacy=False, fmt="tdd")
        logger.info("Completed reload %s at %s (%.4f seconds)", self.id(), time.strftime("%Y %m %d %H:%M:%S", time.localtime()), time.time() - self.__startTime)
        #
        taxId = umP.getMappedId("Q6GZX0", mapName="NCBI-taxon")
        logger.info("TaxId %r", taxId)
        self.assertEqual(taxId, "654924")
        #
        taxId = umP.getMappedIdLegacy("Q6GZX0", mapName="NCBI-taxon")
        logger.info("TaxId %r", taxId)
        # self.assertEqual(taxId, "654924")
        #
        ok = umP.backup(self.__cfgOb, self.__configName)
        self.assertTrue(ok)
        #
        ok = umP.restore(self.__cfgOb, self.__configName)
        self.assertTrue(ok)


def uniProtIdMappingProviderCacheSuite():
    suiteSelect = unittest.TestSuite()
    suiteSelect.addTest(UniProtIdMappingProviderTests("testUniProtIdMappingProviderCachePic"))
    suiteSelect.addTest(UniProtIdMappingProviderTests("testUniProtIdMappingProviderCacheTdd"))
    return suiteSelect


if __name__ == "__main__":
    mySuite = uniProtIdMappingProviderCacheSuite()
    unittest.TextTestRunner(verbosity=2).run(mySuite)
