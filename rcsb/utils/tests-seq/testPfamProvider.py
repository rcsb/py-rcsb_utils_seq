##
# File:    PfamProviderTests.py
# Author:  J. Westbrook
# Date:    18-Feb-2020
#
# Update:
#  26-May-2021 jdw add test for Pfam mapping methods
#
##
"""
Tests utilities to manage access to Pfam data.

"""

__docformat__ = "google en"
__author__ = "John Westbrook"
__email__ = "jwest@rcsb.rutgers.edu"
__license__ = "Apache 2.0"

import logging
import os
import time
import unittest

from rcsb.utils.seq.PfamProvider import PfamProvider

HERE = os.path.abspath(os.path.dirname(__file__))
TOPDIR = os.path.dirname(os.path.dirname(HERE))

logging.basicConfig(level=logging.INFO, format="%(asctime)s [%(levelname)s]-%(module)s.%(funcName)s: %(message)s")
logger = logging.getLogger()


class PfamProviderTests(unittest.TestCase):
    def setUp(self):
        self.__dirPath = os.path.join(os.path.dirname(TOPDIR), "rcsb", "mock-data")
        self.__cachePath = os.path.join(HERE, "test-output")
        self.__startTime = time.time()
        logger.info("Starting %s at %s", self.id(), time.strftime("%Y %m %d %H:%M:%S", time.localtime()))

    def tearDown(self):
        endTime = time.time()
        logger.info("Completed %s at %s (%.4f seconds)\n", self.id(), time.strftime("%Y %m %d %H:%M:%S", time.localtime()), endTime - self.__startTime)

    def testPfamCache(self):
        pP = PfamProvider(cachePath=self.__cachePath, useCache=False)
        ok = pP.testCache()
        self.assertTrue(ok)
        logger.info("Completed %s at %s (%.4f seconds)", self.id(), time.strftime("%Y %m %d %H:%M:%S", time.localtime()), time.time() - self.__startTime)
        # PF18871			HEPN_Toprim_N	HEPN/Toprim N-terminal domain 1
        descr = pP.getDescription("PF18871")
        ok = descr.startswith("HEPN/Toprim N-terminal domain 1")
        self.assertTrue(ok)
        mL = pP.getMapping("1kip")
        logger.debug("mL (%d)", len(mL))
        self.assertGreaterEqual(len(mL), 2)
        vers = pP.getVersion()
        self.assertEqual(vers, "34.0")

    def testPfamCacheFallBack(self):
        pP = PfamProvider(urlTargetPfam="https://rcsb.org/t.txt", urlTargetMapPfam="https://rcsb.org/t.txt", cachePath=self.__cachePath, useCache=False)
        ok = pP.testCache()
        self.assertTrue(ok)
        logger.info("Completed %s at %s (%.4f seconds)", self.id(), time.strftime("%Y %m %d %H:%M:%S", time.localtime()), time.time() - self.__startTime)
        # PF18871			HEPN_Toprim_N	HEPN/Toprim N-terminal domain 1
        descr = pP.getDescription("PF18871")
        ok = descr.startswith("HEPN/Toprim N-terminal domain 1")
        self.assertTrue(ok)
        mL = pP.getMapping("1kip")
        logger.debug("mL (%d)", len(mL))
        self.assertGreaterEqual(len(mL), 2)


def pfamCacheSuite():
    suiteSelect = unittest.TestSuite()
    suiteSelect.addTest(PfamProviderTests("testPfamCache"))
    suiteSelect.addTest(PfamProviderTests("testPfamCacheFallBack"))
    return suiteSelect


if __name__ == "__main__":
    mySuite = pfamCacheSuite()
    unittest.TextTestRunner(verbosity=2).run(mySuite)
