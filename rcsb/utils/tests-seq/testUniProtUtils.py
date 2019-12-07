##
# File:    UniProtUtilsTests.py
# Author:  j. westbrook
# Date:    15-Mar-2019
# Version: 0.001
#
# Update:
##
"""
Test cases for individual and batch fetch of UniProt sequence entries.

    @unittest.skipIf(condition, reason)
    @unittest.skipUnless(condition, reason)
"""

import logging
import os
import unittest

from rcsb.utils.io.MarshalUtil import MarshalUtil
from rcsb.utils.seq.UniProtUtils import UniProtUtils

HERE = os.path.abspath(os.path.dirname(__file__))
TOPDIR = os.path.dirname(os.path.dirname(os.path.dirname(HERE)))

logging.basicConfig(level=logging.INFO, format="%(asctime)s [%(levelname)s]-%(module)s.%(funcName)s: %(message)s")
logger = logging.getLogger()


class UniProtUtilsTests(unittest.TestCase):
    def setUp(self):
        self.__export = False
        self.__mU = MarshalUtil()
        self.__dirPath = os.path.join(os.path.dirname(TOPDIR), "rcsb", "mock-data")
        #
        self.__workPath = os.path.join(HERE, "test-output")
        # Pick up site information from the environment or failover to the development site id.
        #
        self.__unpIdList1 = ["P20937", "P21877", "P22868", "P23832"]
        self.__unpIdList2 = [
            "P29490",
            "P29496",
            "P29498",
            "P29499",
            "P29503",
            "P29506",
            "P29508",
            "P29509",
            "P29525",
            "P29533",
            "P29534",
            "P29547",
            "P29549",
            "P29555",
            "P29557",
            "P29558",
            "P29559",
            "P29563",
            "P29588",
            "P29589",
            "P29590",
            "P29597",
            "P29599",
            "P29600",
            "P29602",
            "P29603",
            "P29617",
            "P29678",
            "P29691",
            "P29715",
            "P29717",
            "P29723",
            "P29724",
            "P29736",
            "P29741",
            "P29745",
            "P29748",
            "P29749",
            "P29752",
            "P29758",
            "P29768",
            "P29803",
            "P29808",
            "P29813",
            "P29827",
            "P29830",
            "P29837",
            "P29838",
            "P29846",
            "P29848",
            "P29882",
            "P29894",
            "P29898",
            "P29899",
            "P29929",
            "P29946",
            "P29957",
            "P29960",
            "P29965",
            "P29966",
            "P29972",
            "P29978",
            "P29986",
            "P29987",
            "P29988",
            "P29989",
            "P29990",
            "P29991",
            "P29994",
        ]
        # self.__unpIdListV = ["P42284-1", "P42284-2", "P42284-3", "P29994-1", "P29994-2", "P29994-3", "P29994-4", "P29994-5", "P29994-6", "P29994-7"]
        self.__unpIdListV = ["P42284-1", "P42284-2", "P42284-3"]

    def testFetchIds(self):
        """ Test individual entry fetch
        """
        try:
            fobj = UniProtUtils(saveText=True)
            for tId in self.__unpIdList1:
                idList = [tId]
                retD, matchD = fobj.fetchList(idList)
                numPrimary, numSecondary, numNone = self.__matchSummary(matchD)
                logger.debug("%d %d %d", numPrimary, numSecondary, numNone)
                #
                rematchD = fobj.rebuildMatchResultIndex(idList, retD)
                self.assertDictEqual(matchD, rematchD)
                #
                self.assertGreaterEqual(len(retD), len(idList))
                if retD and self.__export:
                    fobj.writeUnpXml(os.path.join(self.__workPath, tId + ".xml"))
                    self.__mU.doExport(os.path.join(self.__workPath, tId + ".json"), retD, fmt="json", indent=3)
        except Exception as e:
            logger.exception("Failing with %s", str(e))
            self.fail()

    def testBatchFetch(self):
        """ Test batch entry fetch
        """
        try:
            fobj = UniProtUtils(saveText=False)
            idList = self.__unpIdList1
            retD, matchD = fobj.fetchList(idList)
            numPrimary, numSecondary, numNone = self.__matchSummary(matchD)
            logger.debug("%d %d %d", numPrimary, numSecondary, numNone)
            self.assertGreaterEqual(len(retD), len(idList))
            if retD and self.__export:
                for rId in retD:
                    self.__mU.doExport(os.path.join(self.__workPath, rId + ".json"), retD[rId], fmt="json", indent=3)
        except Exception as e:
            logger.exception("Failing with %s", str(e))
            self.fail()

    def testFetchVariantIds(self):
        """ Test individual variant entry fetch
        """
        try:
            fobj = UniProtUtils(saveText=True)
            for tId in self.__unpIdListV:
                retD, matchD = fobj.fetchList([tId])
                numPrimary, numSecondary, numNone = self.__matchSummary(matchD)
                logger.debug("%d %d %d", numPrimary, numSecondary, numNone)
                self.assertGreaterEqual(len(retD), 1)
                if retD:
                    if self.__export:
                        fobj.writeUnpXml(os.path.join(self.__workPath, tId + ".xml"))
                        self.__mU.doExport(os.path.join(self.__workPath, tId + ".json"), retD, fmt="json", indent=3)
                    #
                    for (eId, eDict) in retD.items():
                        if "db_isoform" in eDict and eId == tId:
                            logger.debug("------ sequence database code  %s has key db_isoform:  %r", eId, eDict["db_isoform"])
                            logger.debug("------ sequence database code  %s sequence length %d", eId, len(eDict["sequence"]))
                            # logger.debug("%s\n", eDict['sequence'])
                        elif eId == tId:
                            logger.debug("------ No matching isoform for %s\n", tId)
                        # for k,v in eDict.items():
                        #    logger.info("%-25s = %s\n", k, v)
                else:
                    logger.info("Fetch failed for id %s\n", tId)
        except Exception as e:
            logger.exception("Failing with %s", str(e))
            self.fail()

    def testBatchFetchVariants(self):
        """  Test batch variant entry fetch
        """
        try:
            fobj = UniProtUtils(saveText=True)
            retD, matchD = fobj.fetchList(self.__unpIdListV)
            numPrimary, numSecondary, numNone = self.__matchSummary(matchD)
            logger.debug("%d %d %d", numPrimary, numSecondary, numNone)
            self.assertGreaterEqual(len(retD), len(self.__unpIdListV))
            if retD and self.__export:
                fobj.writeUnpXml(os.path.join(self.__workPath, "variant-batch-fetch.xml"))
                self.__dumpEntries(retD)
        except Exception as e:
            logger.exception("Failing with %s", str(e))
            self.fail()

    def __dumpEntries(self, retD):
        for (eId, eDict) in retD.items():
            logger.info("------ Entry id %s", eId)
            for k, v in eDict.items():
                logger.info("%-15s = %r", k, v)

    def __matchSummary(self, matchD):
        numPrimary = 0
        numSecondary = 0
        numNone = 0
        for _, mD in matchD.items():
            if mD["matched"] == "primary":
                numPrimary += 1
            elif mD["matched"] == "secondary":
                numSecondary += 1
            else:
                numNone += 1
        logger.info("Matched:  primary:  %d secondary: %d none %d", numPrimary, numSecondary, numNone)
        return numPrimary, numSecondary, numNone


def suiteFetchTests():
    suiteSelect = unittest.TestSuite()
    suiteSelect.addTest(UniProtUtilsTests("testFetchIds"))
    suiteSelect.addTest(UniProtUtilsTests("testBatchFetch"))
    #
    return suiteSelect


def suiteFetchVariantTests():
    suiteSelect = unittest.TestSuite()
    suiteSelect.addTest(UniProtUtilsTests("testFetchVariantIds"))
    suiteSelect.addTest(UniProtUtilsTests("testBatchFetchVariants"))
    #
    return suiteSelect


if __name__ == "__main__":

    #
    mySuite = suiteFetchTests()
    unittest.TextTestRunner(verbosity=2).run(mySuite)
    #
    mySuite = suiteFetchVariantTests()
    unittest.TextTestRunner(verbosity=2).run(mySuite)
#
#
