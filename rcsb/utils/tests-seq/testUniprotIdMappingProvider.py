##
# File:    UniProtIdMappingProviderProviderTests.py
# Author:  J. Westbrook
# Date:    30-Nov-2020
#
# Update:
#
#
##
"""
Tests utilities to manage access to UniProtIdMappingProvider.

"""

__docformat__ = "restructuredtext en"
__author__ = "John Westbrook"
__email__ = "jwest@rcsb.rutgers.edu"
__license__ = "Apache 2.0"

import logging
import os
import platform
import resource
import time
import unittest

from rcsb.utils.seq.UniProtIdMappingProvider import UniProtIdMappingProvider

HERE = os.path.abspath(os.path.dirname(__file__))
TOPDIR = os.path.dirname(os.path.dirname(HERE))

logging.basicConfig(level=logging.INFO, format="%(asctime)s [%(levelname)s]-%(module)s.%(funcName)s: %(message)s")
logger = logging.getLogger()


class UniProtIdMappingProviderTests(unittest.TestCase):
    skipFull = True

    def setUp(self):
        self.__dirPath = os.path.join(os.path.dirname(TOPDIR), "rcsb", "mock-data")
        self.__cachePath = os.path.join(HERE, "test-output")
        self.__startTime = time.time()
        logger.info("Starting %s at %s", self.id(), time.strftime("%Y %m %d %H:%M:%S", time.localtime()))

    def tearDown(self):
        unitS = "MB" if platform.system() == "Darwin" else "GB"
        rusageMax = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
        logger.info("Maximum resident memory size %.4f %s", rusageMax / 10 ** 6, unitS)
        endTime = time.time()
        logger.info("Completed %s at %s (%.4f seconds)", self.id(), time.strftime("%Y %m %d %H:%M:%S", time.localtime()), endTime - self.__startTime)

    @unittest.skipIf(skipFull, "Very long test")
    def testUniProtIdMappingProviderCache(self):
        # umP = UniProtIdMappingProvider(cachePath=self.__cachePath, useCache=True, maxLimit=100)
        umP = UniProtIdMappingProvider(cachePath=self.__cachePath, useCache=True, useLegacy=True)
        logger.info("Completed %s at %s (%.4f seconds)", self.id(), time.strftime("%Y %m %d %H:%M:%S", time.localtime()), time.time() - self.__startTime)
        #
        taxId = umP.getMappedId("Q6GZX0", mapName="NCBI-taxon")
        logger.info("TaxId %r", taxId)
        self.assertEqual(taxId, "654924")
        #
        taxId = umP.getMappedIdLegacy("Q6GZX0", mapName="NCBI-taxon")
        logger.info("TaxId %r", taxId)
        # self.assertEqual(taxId, "654924")


def uniProtIdMappingProviderCacheSuite():
    suiteSelect = unittest.TestSuite()
    suiteSelect.addTest(UniProtIdMappingProviderTests("testUniProtIdMappingProviderCache"))
    return suiteSelect


if __name__ == "__main__":
    mySuite = uniProtIdMappingProviderCacheSuite()
    unittest.TextTestRunner(verbosity=2).run(mySuite)
