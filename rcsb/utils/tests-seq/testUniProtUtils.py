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

    Collaboritive Resarch (CIBR) Proposal by the wwPDB to develop new FACT data services to deliver extreme FAIR data products
            for next generation structural biology to diverse audience of research users
"""

import logging
import os
import unittest

from jsonschema import Draft4Validator
from jsonschema import FormatChecker

from rcsb.utils.io.MarshalUtil import MarshalUtil
from rcsb.utils.seq.UniProtUtils import UniProtUtils

HERE = os.path.abspath(os.path.dirname(__file__))
TOPDIR = os.path.dirname(os.path.dirname(os.path.dirname(HERE)))

logging.basicConfig(level=logging.INFO, format="%(asctime)s [%(levelname)s]-%(module)s.%(funcName)s: %(message)s")
logger = logging.getLogger()


class UniProtUtilsTests(unittest.TestCase):
    def setUp(self):
        self.__export = True
        self.__mU = MarshalUtil()
        self.__dirPath = os.path.join(os.path.dirname(TOPDIR), "rcsb", "mock-data")
        #
        self.__workPath = os.path.join(HERE, "test-output")
        #
        # These settings are for issues with Ubuntu 20.08  (primary site is not reliable)
        self.__usePrimary = True
        self.__retryAltApi = True
        #
        self.__unpIdListMix = ["P20937", "P21877", "P42284-1", "P42284-3", "P29994-2", "P21877"]  # Contains both active and obsoleted entries
        self.__unpIdList1 = ["P20937", "P22868", "P23832", "P21877"]
        self.__unpIdList3 = ["P20937", "P22868", "P23832"]
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
        rawS = """O00422 P00221 Q10423 Q01663 J9VQ51 Q8U004 A0A0H2ZZ96 Q9X1M7 P83876 Q9NPH0 Q74085 Q61466 F1RU54 C7LV29 Q8IHT9 B8XA40 O70439 K0BUH0
                  P23726 Q12748 Q56254 A0A0J9X1Y2 F1LW30 P14625 Q8WV99 Q7YXL2 Q9X250 O07177 P10775 Q38AR9 P48556 A0A0F7RI76 P08107 G9I930 P62399
                  P27254 A8JEA1 E3TFZ1 P13491 P0C249 Q94G94 P78371 Q966D4 Q9KTB7 A0A1L7NZN4 A1SBA1 O15031 G0S0N0 P72223 Q9UUB1 A5U127
                    A0A029IJY9 Q8ZA87 Q93UV0 P20711 A0QWT2 P21453 C5HMM2 Q96DI7 O60760 Q9VVX0 Q7CJ22 Q9REC4 O13516 F2REK9 P0A707 P96724
                    Q8GME2 B3F2X3 A3DK57 O29370 O60494 Q65ZI1 P03886 Q27198 P62375 O27271 A0A0D6J3X3
                    P45628 Q9C471 P96142 Q99720 P70206 Q2RWM6 C7CXJ5 A0QS88 Q83A84 O14793 Q4DSU0 P04806 P48234
                    F8W3X3 P26339 Q9LBQ9 Q9CXY6 P50053 A0A250WER3 A7NNM5 C4B8B9 P08411 P35968 O43291 G4ENZ9 Q01939 C5A086
                    P03436 P30014 Q12314 Q03237 Q08001 P55858 Q05066 Q5ZX93 B1YW99 S5M825 A0A0J9X234 O67413 A0A246CDW5
                    A9PKC6 Q8VVD2 A7RZU9 P06494 Q15059 O96184 Q5FTL8 Q939T0 C9X2N5 W0G557 A0A3F2YM17 P64423 P54573
                    Q926Z8 P9WPY9 A0A059Q5E8 G3X982 B4EEE3 P15428 Q26806 P72391 F6H697 P03925 Q04DI1 A0A095TT41 A0A1S4NYF9 P04075 Q9C5R7
                    O15565 P54652 Q8TAX2 M9MDK9 R4GRT4 A0A0M1WBA0 P11086 Q835Y1 P02752 A7XUK7
                    P49723 Q58801 Q8YY42 Q55670 Q64346 A0A1V3CQ74 Q7X0D9 W5PVD7 E0RXM0 Q5Y7F6 A0A2K3CRG7 P10898 Q57W62 Q03E61 P28804 P29994
                    P16293 O25759 P11236 Q6RFZ0 Q32L40 Q5SLP6 A0A140UGH4 Q85FP8 Q9C3Y6 C3SU37
                    Q8U0N5 D0VWX6 O54050 Q4AE70 Q9H9S0 Q45462 D7PC21 Q5SME6 D0VWR7 P10537 Q8IX01 Q0BRB0 Q3E830 Q58855
                    Q9V1U8 Q8Y5K1 Q52612 P10809 P0AC69 P08905 Q02934 P00378 Q3SYS0 A0A063XHI7 B2XJG5 B7LFT7 Q52SW3
                    Q8BL48 C1DGZ6 G0S4H3 O43602 A0A1S6YJF3 Q6QR64 Q7Z144 F8G0M0 A0R5B5 D0A7Z9 B7N6J5 A6LBR3 P16856
                    P80377 Q7VL95 B8Y5U7 D3KFX5 Q4DIV9 Q92879 O58655 Q84II3 Q9UBU7 P72986 A0A0U4VN94 P32582 P0AC25
                    Q8WUM0 Q9Z8L4 C4ZCS9 P14207 A0A0K8P8E7 P13647 M4MQ92 Q8U0M5 Q9HUK6 Q86FP9 A0A0M3KKX1 Q8TLW1 A0A4V8GZQ8
                    A0A161CFW5 D0VWY5 Q46704 Q8LPB4 Q97Z83 B5BP20 A7ZRI8 Q8GME8 Q00610 A0A0H3JWL8 P38348
                    H9L447 B7MCT6 D1C7H4 E7CH51 Q60I25 Q9I640 G0SHK5 P55776 Q8KRV3 D8IYL4 P70080 Q4QH17 Q56026 P39639
                    Q95W15 P54619 Q9ERI2 P0A185 H9IUR0 B1YQ53 A0QUH3 Q64845 Q9LL85 P09056 K9TLZ5 P35169
                    A5IFX1 P0C0E6 A0A3P3Q1W7 P46436 Q81JF8 A0A482LMF4 A0A0H3JX61 P37362 P52732 Q7T2I5
                    P53974 P00137 Q9R0M6 Q5WFD8 D3E4S5"""
        self.__unpIdListLong = [uId.strip() for uId in rawS.split(" ") if uId]
        self.__jsonSchemaPath = os.path.join(HERE, "test-data", "json-schema-core_uniprot.json")

    def testFetchSequenceList(self):
        """Test fetch UniProt sequence data (FASTA)"""
        try:
            #
            fobj = UniProtUtils(saveText=False)
            # Note: this list contains one obsolete entry
            idList = self.__unpIdList1
            ok, sD = fobj.fetchSequenceList(idList, usePrimary=self.__usePrimary, retryAltApi=self.__retryAltApi)
            logger.debug("sD: %r", sD)
            self.assertFalse(ok)
            self.assertEqual(len(sD), len(idList) - 1)
            if self.__export:
                self.__mU.doExport(os.path.join(self.__workPath, "data-sequences.json"), sD, fmt="json", indent=3)
        except Exception as e:
            logger.exception("Failing with %s", str(e))
            self.fail()

    def testExchangeObject(self):
        """Test fetch exchange objects"""
        try:
            #
            fobj = UniProtUtils(saveText=False)
            idList = self.__unpIdList1
            retD, _ = fobj.fetchList(idList, usePrimary=self.__usePrimary, retryAltApi=self.__retryAltApi)
            exObjD = fobj.reformat(retD, formatType="exchange")
            if exObjD and self.__export:
                for rId in exObjD:
                    self.__mU.doExport(os.path.join(self.__workPath, rId + "-exchange.json"), exObjD[rId], fmt="json", indent=3)
        except Exception as e:
            logger.exception("Failing with %s", str(e))
            self.fail()

    def testValidateExchangeObject(self):
        """Test fetch exchange objects"""
        try:
            #
            sD = self.__mU.doImport(self.__jsonSchemaPath, "json")
            #
            fobj = UniProtUtils(saveText=False)
            idList = self.__unpIdList1
            retD, _ = fobj.fetchList(idList, usePrimary=self.__usePrimary, retryAltApi=self.__retryAltApi)
            #
            exObjD = fobj.reformat(retD, formatType="exchange")
            if exObjD and self.__export:
                for rId in exObjD:
                    self.__mU.doExport(os.path.join(self.__workPath, rId + "-exchange.json"), exObjD[rId], fmt="json", indent=3)
            #
            Draft4Validator.check_schema(sD)
            #
            valInfo = Draft4Validator(sD, format_checker=FormatChecker())
            eCount = 0
            for rId, dD in exObjD.items():
                logger.debug("Uid %s", rId)
                try:
                    cCount = 0
                    for error in sorted(valInfo.iter_errors(dD), key=str):
                        logger.info("%s path %s error: %s", rId, error.path, error.message)
                        logger.debug(">>> failing object is %r", dD)
                        eCount += 1
                        cCount += 1
                    #
                    logger.debug("%s errors count %d", rId, cCount)
                except Exception as e:
                    logger.exception("Validator fails  %s", str(e))
            #
            logger.debug("Total errors count %d", eCount)
        except Exception as e:
            logger.exception("Failing with %s", str(e))
            self.fail()

    def testFetchIds(self):
        """Test individual entry fetch"""
        idList = None
        try:
            fobj = UniProtUtils(saveText=True)
            for tId in self.__unpIdList3:
                idList = [tId]
                retD, matchD = fobj.fetchList(idList, usePrimary=self.__usePrimary, retryAltApi=self.__retryAltApi)
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
            logger.exception("Failing with idList %r %s", idList, str(e))
            self.fail()

    def testBatchFetch(self):
        """Test batch entry fetch"""
        try:
            fobj = UniProtUtils(saveText=False)
            idList = self.__unpIdListMix + self.__unpIdListLong[:100]
            logger.info("idList length %d  unique %d", len(idList), len(set(idList)))
            retD, matchD = fobj.fetchList(idList, usePrimary=self.__usePrimary, retryAltApi=self.__retryAltApi)
            logger.info("IdList %d reference return length %d match length %d", len(idList), len(retD), len(matchD))
            numPrimary, numSecondary, numNone = self.__matchSummary(matchD)
            logger.info("%d %d %d", numPrimary, numSecondary, numNone)
            sumRet = numPrimary + numSecondary + numNone
            logger.info("sumRet returned %d", sumRet)
            self.assertGreaterEqual(sumRet, len(idList) - 1)
            if retD and self.__export:
                for rId in retD:
                    self.__mU.doExport(os.path.join(self.__workPath, rId + ".json"), retD[rId], fmt="json", indent=3)
            #
            logger.info("Test secondary site batch fetch...")
            idListSecondary = idList[0:20]
            retD, matchD = fobj.fetchList(idListSecondary, usePrimary=False, retryAltApi=True, maxChunkSize=10)
            logger.info("idListSecondary %d reference return length %d match length %d", len(idListSecondary), len(retD), len(matchD))
            numPrimary, numSecondary, numNone = self.__matchSummary(matchD)
            logger.debug("%d %d %d", numPrimary, numSecondary, numNone)
            sumRet = numPrimary + numSecondary + numNone
            logger.info("sumRet returned %d", sumRet)
            self.assertGreaterEqual(sumRet, len(idListSecondary) - 1)
            #
        except Exception as e:
            logger.exception("Failing with %s", str(e))
            self.fail()

    def testBatchFetchSecondary(self):
        """Test batch entry fetch (secondary service)"""
        try:
            fobj = UniProtUtils(saveText=False)
            idList = self.__unpIdListLong[:20]
            logger.info("idList length %d  unique %d", len(idList), len(set(idList)))
            try:
                retD, matchD = fobj.fetchList(idList, usePrimary=False, retryAltApi=self.__retryAltApi)
            except Exception as e:
                logger.warning("Fallback secondary service failed with %r", e)
                logger.warning("Retrying with primary service")
                retD, matchD = fobj.fetchList(idList, usePrimary=True, retryAltApi=False)
            logger.info("IdList %d reference return length %d match length %d", len(idList), len(retD), len(matchD))
            numPrimary, numSecondary, numNone = self.__matchSummary(matchD)
            logger.debug("%d %d %d", numPrimary, numSecondary, numNone)
            sumRet = numPrimary + numSecondary + numNone
            logger.info("sumRet returned %d", sumRet)
            self.assertGreaterEqual(sumRet, len(idList) - 1)
            if retD and self.__export:
                for rId in retD:
                    self.__mU.doExport(os.path.join(self.__workPath, rId + ".json"), retD[rId], fmt="json", indent=3)
            #
            retD, matchD = fobj.fetchList(idList, usePrimary=False, retryAltApi=True)
            logger.info("IdList %d reference return length %d match length %d", len(idList), len(retD), len(matchD))
            numPrimary, numSecondary, numNone = self.__matchSummary(matchD)
            logger.debug("%d %d %d", numPrimary, numSecondary, numNone)
            sumRet = numPrimary + numSecondary + numNone
            logger.info("sumRet returned %d", sumRet)
            self.assertGreaterEqual(sumRet, len(idList) - 1)
            #
        except Exception as e:
            logger.exception("Failing with %s", str(e))
            self.fail()

    def testFetchVariantIds(self):
        """Test individual variant entry fetch"""
        try:
            fobj = UniProtUtils(saveText=True)
            for tId in self.__unpIdListV:
                retD, matchD = fobj.fetchList([tId], usePrimary=self.__usePrimary, retryAltApi=self.__retryAltApi)
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
        """Test batch variant entry fetch"""
        try:
            fobj = UniProtUtils(saveText=True)
            retD, matchD = fobj.fetchList(self.__unpIdListV, usePrimary=self.__usePrimary, retryAltApi=self.__retryAltApi)
            numPrimary, numSecondary, numNone = self.__matchSummary(matchD)
            logger.debug("%d %d %d", numPrimary, numSecondary, numNone)
            self.assertGreaterEqual(len(retD), len(self.__unpIdListV))
            if retD and self.__export:
                fobj.writeUnpXml(os.path.join(self.__workPath, "variant-batch-fetch.xml"))
                # self.__dumpEntries(retD)
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

    def testLookup(self):
        """Test lookup gene names"""
        try:
            uUtils = UniProtUtils(saveText=False)
            geneList = ["BCOR"]
            for gene in geneList:
                idList, retCode = uUtils.doLookup([gene], itemKey="Gene_Name")
                logger.info("retCode %r rspList (%d) %r", retCode, len(idList), idList)
                self.assertGreaterEqual(len(idList), 500)
        except Exception as e:
            logger.exception("Failing with %s", str(e))
            self.fail()

    def testGeneLookup(self):
        """Test lookup gene names for human"""
        try:
            uUtils = UniProtUtils(saveText=False)
            geneList = ["BCOR", "BCORL1"]
            for gene in geneList:
                idList, retCode = uUtils.doGeneLookup(gene, 9606)
                logger.info("retCode %r rspList (%d) %r", retCode, len(idList), idList)
                self.assertGreaterEqual(len(idList), 1)
        except Exception as e:
            logger.exception("Failing with %s", str(e))
            self.fail()


def suiteFetchTests():
    suiteSelect = unittest.TestSuite()
    suiteSelect.addTest(UniProtUtilsTests("testLookup"))
    suiteSelect.addTest(UniProtUtilsTests("testFetchIds"))
    suiteSelect.addTest(UniProtUtilsTests("testFetchSequenceList"))
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
