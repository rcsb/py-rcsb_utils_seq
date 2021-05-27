##
#  File:  GlyGenProvider.py
#  Date:  26-May-2021 jdw
#
#  Updates:
#
##
"""
  Fetch glycans and glycoproteins available in the GlyGen.org resource.

"""

import logging
import os.path

from rcsb.utils.io.FileUtil import FileUtil
from rcsb.utils.io.MarshalUtil import MarshalUtil

logger = logging.getLogger(__name__)


class GlyGenProvider(object):
    """Fetch glycans and glycoproteins available in the GlyGen.org resource.

    GlyGen glycan link template -
          https://glygen.org/glycan/G28882EF

    Glycoprotein link template -
          https://www.glygen.org/protein/Q658T7
    """

    def __init__(self, **kwargs):
        #
        self.__cachePath = kwargs.get("cachePath", ".")
        self.__dirPath = os.path.join(self.__cachePath, "glygen")
        useCache = kwargs.get("useCache", True)
        #
        baseUrl = kwargs.get("glygenBasetUrl", "https://data.glygen.org/ln2releases/v-1.8.25/reviewed/")
        fallbackUrl = kwargs.get("glygenFallbackUrl", "https://raw.githubusercontent.com/rcsb/py-rcsb_exdb_assets/master/fall_back/glygen/")
        #
        self.__mU = MarshalUtil(workPath=self.__dirPath)
        self.__glycanD = self.__reloadGlycans(baseUrl, fallbackUrl, self.__dirPath, useCache=useCache)
        self.__glycoproteinD = self.__reloadGlycoproteins(baseUrl, fallbackUrl, self.__dirPath, useCache=useCache)

    def testCache(self, minGlycanCount=20000, minGlycoproteinCount=64000):
        #
        logger.info("GlyGen glycan list (%d) glycoprotein list (%d)", len(self.__glycanD), len(self.__glycoproteinD))
        if self.__glycanD and len(self.__glycanD) > minGlycanCount and self.__glycoproteinD and len(self.__glycoproteinD) > minGlycoproteinCount:
            return True
        return False

    def hasGlycan(self, glyTouCanId):
        try:
            return glyTouCanId in self.__glycanD
        except Exception:
            return False

    def hasGlycoprotein(self, uniProtId):
        try:
            return uniProtId in self.__glycoproteinD
        except Exception:
            return False

    def getGlycans(self):
        return self.__glycanD

    def getGlycoproteins(self):
        return self.__glycoproteinD

    def __reloadGlycans(self, baseUrl, fallbackUrl, dirPath, useCache=True):
        gD = {}
        logger.debug("Using dirPath %r", dirPath)
        self.__mU.mkdir(dirPath)
        #
        myDataPath = os.path.join(dirPath, "glygen-glycan-list.json")
        if useCache and self.__mU.exists(myDataPath):
            gD = self.__mU.doImport(myDataPath, fmt="json")
            logger.debug("GlyGen glycan data length %d", len(gD))
        else:
            logger.debug("Fetch GlyGen glycan data from primary data source %s", baseUrl)
            endPoint = os.path.join(baseUrl, "glycan_masterlist.csv")
            #
            logger.info("Fetch GlyGen glycan data from primary data source %s", endPoint)
            rawPath = os.path.join(dirPath, "glycan_masterlist.csv")
            fU = FileUtil()
            ok = fU.get(endPoint, rawPath)
            logger.debug("Fetch GlyGen glycan data status %r", ok)
            if not ok:
                endPoint = os.path.join(fallbackUrl, "glycan_masterlist.csv")
                ok = fU.get(endPoint, rawPath)
                logger.info("Fetch fallback GlyGen glycan data status %r", ok)
            #
            if ok:
                gD = self.__parseGlycanList(rawPath)
                ok = self.__mU.doExport(myDataPath, gD, fmt="json")
                logger.info("Exported GlyGen glycan list (%d) (%r) %s", len(gD), ok, myDataPath)
            #
        return gD

    def __parseGlycanList(self, filePath):
        gD = {}
        row = None
        try:
            rowL = self.__mU.doImport(filePath, fmt="csv", rowFormat="list")
            logger.debug("Glycan list length (%d)", len(rowL))
            logger.debug("Row 0 %r", rowL[0])
            for row in rowL[1:]:
                gD[row[0]] = row[1]
        except Exception as e:
            logger.exception("Failing for %r (%r) with %s", filePath, row, str(e))
        return gD

    def __reloadGlycoproteins(self, baseUrl, fallbackUrl, dirPath, useCache=True):
        gD = {}
        logger.debug("Using dirPath %r", dirPath)
        self.__mU.mkdir(dirPath)
        #
        myDataPath = os.path.join(dirPath, "glygen-glycoprotein-list.json")
        if useCache and self.__mU.exists(myDataPath):
            gD = self.__mU.doImport(myDataPath, fmt="json")
            logger.debug("GlyGen glycoprotein data length %d", len(gD))
        else:
            for fn in [
                "sarscov1_protein_masterlist.csv",
                "sarscov2_protein_masterlist.csv",
                "hcv1b_protein_masterlist.csv",
                "hcv1a_protein_masterlist.csv",
                "human_protein_masterlist.csv",
                "mouse_protein_masterlist.csv",
                "rat_protein_masterlist.csv",
            ]:
                logger.debug("Fetch GlyGen glycoprotein data from primary data source %s", baseUrl)
                endPoint = os.path.join(baseUrl, fn)
                #
                logger.debug("Fetch GlyGen glycoprotein data from primary data source %s", endPoint)
                rawPath = os.path.join(dirPath, fn)
                fU = FileUtil()
                ok = fU.get(endPoint, rawPath)
                logger.debug("Fetch GlyGen glycoprotein data status %r", ok)
                if not ok:
                    endPoint = os.path.join(fallbackUrl, fn)
                    ok = fU.get(endPoint, rawPath)
                    logger.info("Fetch fallback GlyGen data status %r", ok)
                #
                if ok:
                    tD = self.__parseGlycoproteinList(rawPath)
                    gD.update(tD)
            #
            ok = self.__mU.doExport(myDataPath, gD, fmt="json")
            logger.info("Exported GlyGen glycoprotein list (%d) (%r) %s", len(gD), ok, myDataPath)
        #
        return gD

    def __parseGlycoproteinList(self, filePath):
        gD = {}
        try:
            rowL = self.__mU.doImport(filePath, fmt="csv", rowFormat="list")
            for row in rowL[1:]:
                ff = row[0].split("-")
                gD[ff[0]] = ff[1]
        except Exception as e:
            logger.exception("Failing for %r with %s", filePath, str(e))
        return gD
