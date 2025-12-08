##
# File:    GlyGenProviderTests.py
# Date:    26-May-2021  JDW
#
# Updates:
#
##
"""
Test cases for operations that fetch glycans and glycoproteins from GlyGen.org -

"""

import logging
import os
import time
import unittest

from importlib.metadata import version as get_package_version
from rcsb.utils.seq.GlyGenProvider import GlyGenProvider

__version__ = get_package_version("rcsb.utils.seq")

HERE = os.path.abspath(os.path.dirname(__file__))
TOPDIR = os.path.dirname(os.path.dirname(HERE))

logging.basicConfig(level=logging.INFO, format="%(asctime)s [%(levelname)s]-%(module)s.%(funcName)s: %(message)s")
logger = logging.getLogger()


class GlyGenProviderTests(unittest.TestCase):
    def setUp(self):
        self.__workPath = os.path.join(HERE, "test-output")
        self.__dataPath = os.path.join(HERE, "test-data")
        self.__cachePath = os.path.join(self.__workPath, "CACHE")
        #
        self.__startTime = time.time()
        logger.debug("Running tests on version %s", __version__)
        logger.debug("Starting %s at %s", self.id(), time.strftime("%Y %m %d %H:%M:%S", time.localtime()))

    def tearDown(self):
        endTime = time.time()
        logger.debug("Completed %s at %s (%.4f seconds)", self.id(), time.strftime("%Y %m %d %H:%M:%S", time.localtime()), endTime - self.__startTime)

    def testGetGlyGenData(self):
        """Load GlyGen data sets"""
        try:
            ggP = GlyGenProvider(cachePath=self.__cachePath, useCache=False)
            ok = ggP.testCache()
            self.assertTrue(ok)
            gD = ggP.getGlycans()
            self.assertGreaterEqual(len(gD), 30000)
            for gId in gD:
                ok = ggP.hasGlycan(gId)
                self.assertTrue(ok)
            gD = ggP.getGlycoproteins()
            self.assertGreaterEqual(len(gD), 64000)
            #
            for gId in gD:
                ok = ggP.hasGlycoprotein(gId)
                self.assertTrue(ok)

        except Exception as e:
            logger.exception("Failing with %s", str(e))
            self.fail()

    def testGetGlyGenDataFallback(self):
        """Load GlyGen data sets using fallback site"""
        try:
            ggP = GlyGenProvider(glygenBasetUrl="https://rcsb.org/t", cachePath=self.__cachePath, useCache=False)
            ok = ggP.testCache()
            self.assertTrue(ok)
            gD = ggP.getGlycans()
            self.assertGreaterEqual(len(gD), 30000)
            for gId in gD:
                ok = ggP.hasGlycan(gId)
                self.assertTrue(ok)
            gD = ggP.getGlycoproteins()
            self.assertGreaterEqual(len(gD), 64000)
            #
            for gId in gD:
                ok = ggP.hasGlycoprotein(gId)
                self.assertTrue(ok)

        except Exception as e:
            logger.exception("Failing with %s", str(e))
            self.fail()


def readGlyGenData():
    suiteSelect = unittest.TestSuite()
    suiteSelect.addTest(GlyGenProviderTests("testGetGlyGenData"))
    suiteSelect.addTest(GlyGenProviderTests("testGetGlyGenDataFallback"))
    return suiteSelect


if __name__ == "__main__":
    mySuite = readGlyGenData()
    unittest.TextTestRunner(verbosity=2).run(mySuite)
