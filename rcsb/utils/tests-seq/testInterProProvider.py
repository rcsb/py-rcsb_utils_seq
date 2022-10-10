##
# File:    InterProProviderTests.py
# Author:  J. Westbrook
# Date:    18-Feb-2020
#
# Update:
#
#
##
"""
Tests utilities to manage access to InterPro data.

"""

__docformat__ = "google en"
__author__ = "John Westbrook"
__email__ = "jwest@rcsb.rutgers.edu"
__license__ = "Apache 2.0"

import logging
import os
import time
import unittest

from rcsb.utils.seq.InterProProvider import InterProProvider

HERE = os.path.abspath(os.path.dirname(__file__))
TOPDIR = os.path.dirname(os.path.dirname(HERE))

logging.basicConfig(level=logging.INFO, format="%(asctime)s [%(levelname)s]-%(module)s.%(funcName)s: %(message)s")
logger = logging.getLogger()


class InterProProviderTests(unittest.TestCase):
    def setUp(self):
        self.__dirPath = os.path.join(os.path.dirname(TOPDIR), "rcsb", "mock-data")
        self.__cachePath = os.path.join(HERE, "test-output")
        self.__startTime = time.time()
        logger.info("Starting %s at %s", self.id(), time.strftime("%Y %m %d %H:%M:%S", time.localtime()))

    def tearDown(self):
        endTime = time.time()
        logger.info("Completed %s at %s (%.4f seconds)\n", self.id(), time.strftime("%Y %m %d %H:%M:%S", time.localtime()), endTime - self.__startTime)

    def testInterProCache(self):
        ipP = InterProProvider(cachePath=self.__cachePath, useCache=False, useFallBack=False)
        logger.info("Completed %s at %s (%.4f seconds)", self.id(), time.strftime("%Y %m %d %H:%M:%S", time.localtime()), time.time() - self.__startTime)
        # IPR041653	Repeat	Importin repeat 4
        idCode = "IPR041653"
        txt = ipP.getDescription(idCode)
        ok = txt.startswith("Importin repeat 4")
        logger.info("InterPro %s description %s", idCode, txt)
        self.assertTrue(ok)
        #
        txt = ipP.getType(idCode)
        logger.info("InterPro %s type %s", idCode, txt)
        ok = txt.startswith("Repeat")
        self.assertTrue(ok)

        idCode = "IPR023678"
        linL = ipP.getLineage(idCode)
        logger.debug("lin %r", linL)
        self.assertEqual(len(linL), 4)
        tL = ipP.getTreeNodeList()
        logger.debug("tree list %r", tL[:50])
        logger.info("tree node list length %d", len(tL))
        self.assertGreaterEqual(len(tL), 30000)

    def testInterProCacheFallBack(self):
        """Test case for utilizing fallback data when provider fails to retreive from given urlTarget.
        (You should expect to see multiple errors here, prior to resorting to fall back.)
        """
        ipP = InterProProvider(urlTargetInterPro="https://rcsb.org/t.txt", cachePath=self.__cachePath, useCache=False)
        logger.info("Completed %s at %s (%.4f seconds)", self.id(), time.strftime("%Y %m %d %H:%M:%S", time.localtime()), time.time() - self.__startTime)
        # IPR041653	Repeat	Importin repeat 4
        idCode = "IPR041653"
        txt = ipP.getDescription(idCode)
        ok = txt.startswith("Importin repeat 4")
        logger.info("InterPro %s description %s", idCode, txt)
        self.assertTrue(ok)
        #
        txt = ipP.getType(idCode)
        logger.info("InterPro %s type %s", idCode, txt)
        ok = txt.startswith("Repeat")
        self.assertTrue(ok)
        #
        idCode = "IPR023678"
        linL = ipP.getLineage(idCode)
        logger.debug("lin %r", linL)
        self.assertEqual(len(linL), 4)


def interProCacheSuite():
    suiteSelect = unittest.TestSuite()
    suiteSelect.addTest(InterProProviderTests("testInterProCache"))
    suiteSelect.addTest(InterProProviderTests("testInterProCacheFallBack"))
    return suiteSelect


if __name__ == "__main__":
    mySuite = interProCacheSuite()
    unittest.TextTestRunner(verbosity=2).run(mySuite)
