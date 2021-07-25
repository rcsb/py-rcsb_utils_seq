##
# File:    GlycanProviderTests.py
# Author:  J. Westbrook
# Date:    28-May-2021
#
# Update:
#
#
##
"""
Tests for accessors for managing Glycan extracted annotations. (read-only version)

"""

__docformat__ = "google en"
__author__ = "John Westbrook"
__email__ = "jwest@rcsb.rutgers.edu"
__license__ = "Apache 2.0"

import logging
import os
import unittest

from rcsb.utils.io.FileUtil import FileUtil
from rcsb.utils.seq.GlycanProvider import GlycanProvider
from rcsb.utils.config.ConfigUtil import ConfigUtil

HERE = os.path.abspath(os.path.dirname(__file__))
TOPDIR = os.path.dirname(os.path.dirname(os.path.dirname(HERE)))

logging.basicConfig(level=logging.INFO, format="%(asctime)s [%(levelname)s]-%(module)s.%(funcName)s: %(message)s")
logger = logging.getLogger()


class GlycanProviderTests(unittest.TestCase):
    def setUp(self):
        self.__cachePath = os.path.join(HERE, "test-output", "CACHE")
        self.__dataPath = os.path.join(HERE, "test-data")
        #
        self.__mockTopPath = os.path.join(TOPDIR, "rcsb", "mock-data")
        self.__configPath = os.path.join(TOPDIR, "rcsb", "mock-data", "config", "dbload-setup-example.yml")

    def tearDown(self):
        pass

    def testGlycanMapping(self):

        minCount = 12
        configName = "site_info_configuration"
        cfgOb = ConfigUtil(configPath=self.__configPath, defaultSectionName=configName, mockTopPath=self.__mockTopPath)
        #
        # mock the stashed data for this provider.
        basePath = cfgOb.get("_STASH_SERVER_BASE_PATH", sectionName=configName)
        fU = FileUtil()
        fn = "glycan.tar.gz"
        bundlePath = os.path.join(self.__dataPath, fn)
        stashBundlePath = os.path.join(self.__cachePath, basePath, fn)
        fU.put(bundlePath, stashBundlePath)
        #
        gP = GlycanProvider(cachePath=self.__cachePath, useCache=False)
        ok = gP.testCache(minCount=0)
        self.assertTrue(ok)
        #
        ok = gP.restore(cfgOb, configName)
        self.assertTrue(ok)
        #
        ok = gP.reload()
        self.assertTrue(ok)
        ok = gP.testCache(minCount=minCount)
        self.assertTrue(ok)
        #
        riD = gP.getIdentifiers()
        logger.info("riD (%d)", len(riD))
        self.assertTrue(ok)
        #
        ok = gP.restore(cfgOb, configName)
        self.assertTrue(ok)
        #
        ok = gP.reload()
        self.assertTrue(ok)
        ok = gP.testCache(minCount=minCount)
        self.assertTrue(ok)

    @unittest.skip("Internal test")
    def testGlycanStashRemote(self):
        minCount = 12
        configName = "site_info_remote_configuration"
        cfgOb = ConfigUtil(configPath=self.__configPath, defaultSectionName=configName, mockTopPath=self.__mockTopPath)
        #
        gP = GlycanProvider(cachePath=self.__cachePath, useCache=False)
        ok = gP.testCache(minCount=0)
        self.assertTrue(ok)
        #
        ok = gP.restore(cfgOb, configName)
        self.assertTrue(ok)
        #
        ok = gP.reload()
        self.assertTrue(ok)
        ok = gP.testCache(minCount=minCount)
        self.assertTrue(ok)
        #
        riD = gP.getIdentifiers()
        logger.info("riD (%d)", len(riD))
        self.assertTrue(ok)
        #
        ok = gP.restore(cfgOb, configName)
        self.assertTrue(ok)
        #
        ok = gP.reload()
        self.assertTrue(ok)
        ok = gP.testCache(minCount=minCount)
        self.assertTrue(ok)
        #
        riD = gP.getIdentifiers()
        logger.info("riD (%d)", len(riD))
        self.assertTrue(ok)
        #


def glycanMappingSuite():
    suiteSelect = unittest.TestSuite()
    suiteSelect.addTest(GlycanProviderTests("testGlycanMapping"))
    return suiteSelect


if __name__ == "__main__":

    mySuite = glycanMappingSuite()
    unittest.TextTestRunner(verbosity=2).run(mySuite)
