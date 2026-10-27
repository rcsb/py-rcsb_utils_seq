##
# File:    GlycanProviderTests.py
# Author:  J. Westbrook
# Date:    28-May-2021
#
# Update:
#  20-Oct-2026 dwp Add test for getGlycanIdentifier() with short and extended entry IDs
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
from rcsb.utils.io.MarshalUtil import MarshalUtil
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

    def testGetGlycanIdentifier(self):
        """Test glycan identifier lookups for short and extended entry IDs (keys as stored in the identifier map)."""
        cachePath = os.path.join(HERE, "test-output", "CACHE-glycan-ids")
        idD = {
            "1ABC_2": {"glyTouCanId": "G00001AA", "prdId": "PRD_000001"},
            "pdb_00001abc_3": {"glyTouCanId": "G00002BB", "prdId": None},
            "pdb_10000xyz_1": {"glyTouCanId": "G00003CC", "prdId": None},
        }
        mU = MarshalUtil()
        ok = mU.doExport(os.path.join(cachePath, "glycan", "branched_entity_glycan_identifier_map.json"), {"version": "0.50", "created": "", "identifiers": idD}, fmt="json")
        self.assertTrue(ok)
        gP = GlycanProvider(cachePath=cachePath, useCache=True)
        self.assertTrue(gP.testCache(minCount=3))
        # Short and extended IDs, including an extended ID with no short form
        self.assertEqual(gP.getGlycanIdentifier("1ABC", "2", "glyTouCanId"), "G00001AA")
        self.assertEqual(gP.getGlycanIdentifier("1ABC", "2", "prdId"), "PRD_000001")
        self.assertEqual(gP.getGlycanIdentifier("pdb_00001abc", "3", "glyTouCanId"), "G00002BB")
        self.assertEqual(gP.getGlycanIdentifier("pdb_10000xyz", "1", "glyTouCanId"), "G00003CC")
        # Misses: unknown entity, ID type not present, or a key in a different format/case than stored
        self.assertIsNone(gP.getGlycanIdentifier("1ABC", "9", "glyTouCanId"))
        self.assertIsNone(gP.getGlycanIdentifier("pdb_00001abc", "3", "prdId"))
        self.assertIsNone(gP.getGlycanIdentifier("pdb_00001abc", "2", "glyTouCanId"))
        self.assertIsNone(gP.getGlycanIdentifier("1abc", "2", "glyTouCanId"))

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
    suiteSelect.addTest(GlycanProviderTests("testGetGlycanIdentifier"))
    return suiteSelect


if __name__ == "__main__":

    mySuite = glycanMappingSuite()
    unittest.TextTestRunner(verbosity=2).run(mySuite)
