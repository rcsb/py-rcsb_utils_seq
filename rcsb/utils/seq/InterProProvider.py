##
# File:    InterProProvider.py
# Date:    18-Feb-2020
#
##

import logging
import os
import sys

from rcsb.utils.io.FileUtil import FileUtil
from rcsb.utils.io.MarshalUtil import MarshalUtil

logger = logging.getLogger(__name__)


class InterProProvider(object):
    """ Manage mappings of InterPro identifiers to description and parent/child relationships

    """

    def __init__(self, **kwargs):
        urlTargetInterPro = kwargs.get("urlTargetInterPro", "ftp://ftp.ebi.ac.uk/pub/databases/interpro/current/entry.list")
        urlTargetInterProFB = "https://github.com/rcsb/py-rcsb_exdb_assets/raw/master/fall_back/InterPro/entry.list"
        urlTargetInterProParent = kwargs.get("urlTargetInterPro", "ftp://ftp.ebi.ac.uk/pub/databases/interpro/current/ParentChildTreeFile.txt")
        urlTargetInterProParentFB = "https://github.com/rcsb/py-rcsb_exdb_assets/raw/master/fall_back/InterPro/ParentChildTreeFile.txt"
        cachePath = kwargs.get("cachePath", ".")
        dirPath = os.path.join(cachePath, "interPro")
        useCache = kwargs.get("useCache", True)
        #
        self.__mU = MarshalUtil(workPath=dirPath)
        self.__interProD, self.__interProParentD = self.__rebuildCache(
            urlTargetInterPro, urlTargetInterProFB, urlTargetInterProParent, urlTargetInterProParentFB, dirPath, useCache
        )

    def getDescription(self, interProId):
        ret = None
        try:
            ret = self.__interProD[interProId]["description"]
        except Exception:
            pass
        return ret

    def getType(self, interProId):
        ret = None
        try:
            ret = self.__interProD[interProId]["type"]
        except Exception:
            pass
        return ret

    def testCache(self):
        # Check length ...
        logger.info("Length InterPro %d", len(self.__interProD))
        return len(self.__interProD) > 1000

    #
    def __rebuildCache(self, urlTargetInterPro, urlTargetInterProFB, urlTargetInterProParent, urlTargetInterProParentFB, dirPath, useCache):
        fmt = "json"
        ext = fmt if fmt == "json" else "pic"
        interProDataPath = os.path.join(dirPath, "interPro-data.%s" % ext)
        #
        logger.debug("Using cache data path %s", dirPath)
        self.__mU.mkdir(dirPath)
        #
        if useCache and self.__mU.exists(interProDataPath):
            rD = self.__mU.doImport(interProDataPath, fmt=fmt)
            interProD = rD["index"]
            interProParentD = rD["parents"]
            logger.debug("InterPro index length %d parent length %d", len(interProD), len(interProParentD))
        else:
            # ------
            fU = FileUtil()
            logger.info("Fetch data from source %s in %s", urlTargetInterPro, dirPath)
            fp = os.path.join(dirPath, fU.getFileName(urlTargetInterPro))
            ok = fU.get(urlTargetInterPro, fp)
            if not ok:
                fp = os.path.join(dirPath, fU.getFileName(urlTargetInterProFB))
                ok = fU.get(urlTargetInterProFB, fp)
                logger.info("Fetch data fallback fetch status is %r", ok)
            interProD = self.__getInterProIndex(fp)

            logger.info("Caching %d in %s status %r", len(interProD), interProDataPath, ok)
            # ------
            logger.info("Fetch data from source %s in %s", urlTargetInterPro, dirPath)
            fp = os.path.join(dirPath, fU.getFileName(urlTargetInterProParent))
            ok = fU.get(urlTargetInterProParent, fp)
            if not ok:
                fp = os.path.join(dirPath, fU.getFileName(urlTargetInterProParentFB))
                ok = fU.get(urlTargetInterProParentFB, fp)
                logger.info("Fetch data fallback fetch status is %r", ok)
            interProParentD = self.__getInterProParents(fp)
            #
            ok = self.__mU.doExport(interProDataPath, {"index": interProD, "parents": interProParentD}, fmt=fmt)
        #
        return interProD, interProParentD

    def getLineage(self, idCode):
        pList = []
        try:
            pList.append(idCode)
            pt = self.getParentId(idCode)
            while (pt is not None) and (pt != 1):
                pList.append(pt)
                pt = self.getParentId(pt)
        except Exception as e:
            logger.exception("Failing with %s", str(e))
        #
        pList.reverse()
        return pList

    def getLineageWithNames(self, idCode):
        linL = []
        try:
            idCodeL = self.getLineage(idCode)
            for ii, idCode in enumerate(idCodeL, 1):
                linL.append((idCode, self.getDescription(idCode), ii))
        except Exception as e:
            logger.exception("Failing with %s", str(e))
        return linL

    def getParentId(self, idCode):
        try:
            return self.__interProParentD[idCode]
        except Exception:
            pass
        return None

    def getTreeNodeList(self, filterD=None):
        dL = []
        try:
            for idCode, _ in self.__interProD.items():
                if filterD and idCode not in filterD:
                    continue
                displayName = self.getDescription(idCode)
                pId = self.getParentId(idCode)
                linL = self.getLineage(idCode)
                #
                if pId is None:
                    dD = {"id": idCode, "name": displayName, "depth": 0}
                else:
                    dD = {"id": idCode, "name": displayName, "parents": [pId], "depth": len(linL) - 1}
                dL.append(dD)
        except Exception as e:
            logger.exception("Failing with %s", str(e))
        return dL

    def __getInterProParents(self, filePath):
        """Read the InterPro parent hierarchy and return a dictionary parent ids.

        Args:
            filePath (str): path to InterPro parent/child hierachy

        Returns:
            dict: {idCode: parentIdCode or None}
        """
        interProParentD = {}
        lineL = self.__mU.doImport(filePath, fmt="list")
        stack = []
        for line in lineL:
            content = line.rstrip()  # drop \n
            row = content.split("--")
            ff = row[-1].split("::")
            tS = ff[0].strip()
            # stack[:] = stack[: len(row) - 1] + [row[-1]]
            stack[:] = stack[: len(row) - 1] + [tS]
            for ii, idCode in enumerate(stack):
                interProParentD[idCode] = None if ii == 0 else stack[ii - 1]
            logger.debug("Lineage %r", "\t".join(stack))
        #
        return interProParentD

    def __getInterProIndex(self, filePath):
        """Read CSV file of InterPro accessions and descriptions

        Args:
            filePath (str): path to InterPro accession/description csv file

        Returns:
            dict: {idCode: description}
        """

        interProD = {}
        encodingD = {"encoding": "ascii"} if sys.version_info[0] < 3 else {}
        rowL = self.__mU.doImport(filePath, fmt="tdd", rowFormat="list", **encodingD)
        for row in rowL:
            try:
                interProId = row[0].strip().upper()
                interProType = row[1].strip()
                descr = row[2].strip()
                interProD[interProId] = {"description": descr, "type": interProType}
            except Exception:
                pass
        #
        return interProD
