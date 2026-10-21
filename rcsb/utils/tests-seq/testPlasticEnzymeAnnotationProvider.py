##
# File:    testPlasticEnzymeAnnotationProvider.py
# Date:    21-Oct-2026
# Version: 0.001
#
# Update:
#
##
"""
Tests for utilities managing plastic-degrading enzyme annotation data (PAZy and PlasticDB).

The tests run against the downloaded source files in this directory (pazy_proteins.json, degraders_list.tsv,
pazy_proteins.fasta and PlasticDB.fasta), and check the results against the contents of those files and of
the generated select-data files (rather than against hardcoded values).
"""

__docformat__ = "google en"
__license__ = "Apache 2.0"

import logging
import os
import platform
import resource
import shutil
import time
import unittest

from rcsb.utils.io.MarshalUtil import MarshalUtil
from rcsb.utils.seq.PlasticEnzymeAnnotationProvider import PlasticEnzymeAnnotationProvider

HERE = os.path.abspath(os.path.dirname(__file__))

logging.basicConfig(level=logging.INFO)
logger = logging.getLogger()


class PlasticEnzymeAnnotationProviderTests(unittest.TestCase):
    def setUp(self):
        self.__cachePath = os.path.join(HERE, "test-output", "CACHE")
        self.__dirPath = os.path.join(self.__cachePath, "plastic-enzyme-annotations")
        # Start each test with an empty cache
        shutil.rmtree(self.__cachePath, ignore_errors=True)
        #
        self.__pazyJsonPath = os.path.join(HERE, "pazy_proteins.json")
        self.__plasticDBTsvPath = os.path.join(HERE, "degraders_list.tsv")
        self.__pazyFastaPath = os.path.join(HERE, "pazy_proteins.fasta")
        self.__plasticDBFastaPath = os.path.join(HERE, "PlasticDB.fasta")
        for filePath in [self.__pazyJsonPath, self.__plasticDBTsvPath, self.__pazyFastaPath, self.__plasticDBFastaPath]:
            if not os.path.exists(filePath):
                self.skipTest("Source data file not present: %s" % filePath)
        #
        self.__mU = MarshalUtil()
        self.__startTime = time.time()
        logger.info("Starting %s at %s", self.id(), time.strftime("%Y %m %d %H:%M:%S", time.localtime()))

    def tearDown(self):
        unitS = "MB" if platform.system() == "Darwin" else "GB"
        rusageMax = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
        logger.info("Maximum resident memory size %.4f %s", rusageMax / 10**6, unitS)
        endTime = time.time()
        logger.info("Completed %s at %s (%.4f seconds)", self.id(), time.strftime("%Y %m %d %H:%M:%S", time.localtime()), endTime - self.__startTime)

    # --- helpers
    def __getProvider(self, **kwargs):
        # The source files in this directory are used in place of downloading from PAZy and PlasticDB
        return PlasticEnzymeAnnotationProvider(cachePath=self.__cachePath, pazyUrl=self.__pazyJsonPath, plasticDBUrl=self.__plasticDBTsvPath, **kwargs)

    def __getDataPath(self, fileName):
        return os.path.join(self.__dirPath, fileName)

    def __getSelectData(self, fileName):
        """Return the "data" section of a generated select-data file."""
        fD = self.__mU.doImport(self.__getDataPath(fileName), fmt="json")
        self.assertIn("version", fD)
        self.assertIn("created", fD)
        self.assertIn("data", fD)
        return fD["data"]

    def __getFastaRecords(self, fileName):
        """Return a list of (header, sequence) tuples from a generated fasta file."""
        recordL = []
        with open(self.__getDataPath(fileName), "r", encoding="utf-8") as ifh:
            for line in ifh.read().splitlines():
                if line.startswith(">"):
                    recordL.append([line[1:], ""])
                else:
                    recordL[-1][1] += line
        return [tuple(r) for r in recordL]

    def __checkFastaMatchesSelectData(self, fastaFileName, pazyD, plasticDBD):
        """The fasta has one record per select-data entry (PAZy first, then PlasticDB), with header >{id}|enzymeId|{sourceDb}|sourceDb"""
        expectedL = []
        for dataD in (pazyD, plasticDBD):
            for enzymeId, eD in dataD.items():
                expectedL.append(("%s|enzymeId|%s|sourceDb" % (enzymeId, eD["source_db"]), eD["sequence"]))
        self.assertEqual(self.__getFastaRecords(fastaFileName), expectedL)

    # --- tests
    def testFetchPlasticEnzymeData(self):
        pzP = self.__getProvider(useCache=False)
        ok = pzP.testCache()
        self.assertFalse(ok)
        ok = pzP.fetchPlasticEnzymeData()
        self.assertTrue(ok)
        #
        for fileName, sourceDb in [("pazy-select-data.json", "PAZy"), ("plasticdb-select-data.json", "PlasticDB")]:
            dataD = self.__getSelectData(fileName)
            self.assertGreater(len(dataD), 100)
            for enzymeId, eD in dataD.items():
                self.assertEqual(eD["source_db"], sourceDb)
                self.assertTrue(eD["sequence"], msg="%s %s" % (sourceDb, enzymeId))
                self.assertTrue(eD["name"], msg="%s %s" % (sourceDb, enzymeId))
                self.assertIsInstance(eD["substrates"], list)

    def testExportFasta(self):
        pzP = self.__getProvider(useCache=False)
        # Fasta cannot be exported before the data has been fetched
        ok = pzP.exportFasta()
        self.assertFalse(ok)
        ok = pzP.fetchPlasticEnzymeData()
        self.assertTrue(ok)
        ok = pzP.exportFasta()
        self.assertTrue(ok)
        #
        self.__checkFastaMatchesSelectData("plastic-enzyme-sequences.fa", self.__getSelectData("pazy-select-data.json"), self.__getSelectData("plasticdb-select-data.json"))

    def testExportFastaFromFastaFiles(self):
        pzP = self.__getProvider(useCache=False)
        ok = pzP.exportFastaFromFastaFiles(self.__pazyFastaPath, self.__plasticDBFastaPath)
        self.assertTrue(ok)
        #
        pazyD = self.__getSelectData("pazy-fasta-select-data.json")
        plasticDBD = self.__getSelectData("plasticdb-fasta-select-data.json")
        self.assertGreater(len(pazyD), 100)
        self.assertGreater(len(plasticDBD), 100)
        self.__checkFastaMatchesSelectData("fa-plastic-enzymes-sequences.fa", pazyD, plasticDBD)
        #
        # Both workflows give the same set of sequences
        ok = pzP.fetchPlasticEnzymeData()
        self.assertTrue(ok)
        ok = pzP.exportFasta()
        self.assertTrue(ok)
        self.assertEqual(dict(self.__getFastaRecords("fa-plastic-enzymes-sequences.fa")), dict(self.__getFastaRecords("plastic-enzyme-sequences.fa")))

    def testReloadAndTestCache(self):
        pzP = self.__getProvider(useCache=False)
        ok = pzP.testCache()
        self.assertFalse(ok)
        ok = pzP.fetchPlasticEnzymeData()
        self.assertTrue(ok)
        # The select-data files are not loaded until reloaded; useCache=False loads nothing
        pzP = self.__getProvider(useCache=False)
        self.assertFalse(pzP.testCache())
        #
        # The annotation data file (generated by buildAnnotationList) is not yet available, so the cache is not complete
        pzP = self.__getProvider(useCache=True)
        self.assertFalse(pzP.testCache())
        #
        # Create an annotation data file (one entry per enzyme) in place of the output of buildAnnotationList
        aD = {}
        for fileName in ["pazy-select-data.json", "plasticdb-select-data.json"]:
            for enzymeId, eD in self.__getSelectData(fileName).items():
                aD["%s_1" % enzymeId] = {"entry_id": enzymeId, "entity_id": "1", "annotation_id": enzymeId, "source_db": eD["source_db"]}
        self.assertGreater(len(aD), 100)
        ok = self.__mU.doExport(self.__getDataPath("plastic-enzyme-annotation-data.json"), {"version": "1.0", "created": "test", "annotations": aD}, fmt="json", indent=3)
        self.assertTrue(ok)
        #
        ok = pzP.reload()
        self.assertTrue(ok)
        ok = pzP.testCache()
        self.assertTrue(ok)
        pzP = self.__getProvider(useCache=True)
        ok = pzP.testCache()
        self.assertTrue(ok)
        # useCache=False still loads nothing
        pzP = self.__getProvider(useCache=False)
        ok = pzP.testCache()
        self.assertFalse(ok)


def buildPlasticEnzymeAnnotationTests():
    suiteSelect = unittest.TestSuite()
    suiteSelect.addTest(PlasticEnzymeAnnotationProviderTests("testFetchPlasticEnzymeData"))
    suiteSelect.addTest(PlasticEnzymeAnnotationProviderTests("testExportFasta"))
    suiteSelect.addTest(PlasticEnzymeAnnotationProviderTests("testExportFastaFromFastaFiles"))
    suiteSelect.addTest(PlasticEnzymeAnnotationProviderTests("testReloadAndTestCache"))
    return suiteSelect


if __name__ == "__main__":
    mySuite = buildPlasticEnzymeAnnotationTests()
    unittest.TextTestRunner(verbosity=2).run(mySuite)
