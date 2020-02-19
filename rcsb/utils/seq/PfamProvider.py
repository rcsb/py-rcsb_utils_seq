##
# File:    PfamProvider.py
# Date:    18-Feb-2020
#
##

import logging
import os
import sys

from rcsb.utils.io.FileUtil import FileUtil
from rcsb.utils.io.MarshalUtil import MarshalUtil

logger = logging.getLogger(__name__)


class PfamProvider(object):
    """ Manage an index of Pfam identifier to description mappings.
    """

    def __init__(self, **kwargs):
        urlTargetPfam = kwargs.get("urlTargetPfam", "ftp://ftp.ebi.ac.uk/pub/databases/Pfam/current_release/Pfam-A.clans.tsv.gz")
        urlTargetPfamFB = "https://github.com/rcsb/py-rcsb_exdb_assets/raw/master/fall_back/Pfam/Pfam-A.clans.tsv.gz"
        cachePath = kwargs.get("cachePath", ".")
        dirPath = os.path.join(cachePath, "pfam")
        useCache = kwargs.get("useCache", True)
        #
        self.__mU = MarshalUtil(workPath=dirPath)
        self.__pfamD = self.__rebuildCache(urlTargetPfam, urlTargetPfamFB, dirPath, useCache)

    def getDescription(self, pfamId):
        descr = None
        try:
            descr = self.__pfamD[pfamId]
        except Exception:
            pass
        return descr

    def testCache(self):
        # Check length ...
        logger.info("Length Pfam %d", len(self.__pfamD))
        return len(self.__pfamD) > 1000

    #
    def __rebuildCache(self, urlTargetPfam, urlTargetPfamFB, dirPath, useCache):
        fmt = "json"
        ext = fmt if fmt == "json" else "pic"
        pfamDataPath = os.path.join(dirPath, "pfam-data.%s" % ext)
        #
        logger.debug("Using cache data path %s", dirPath)
        self.__mU.mkdir(dirPath)
        #
        if useCache and self.__mU.exists(pfamDataPath):
            pfamD = self.__mU.doImport(pfamDataPath, fmt=fmt)
            logger.debug("Pfam data length %d", len(pfamD))
        else:
            # ------
            fU = FileUtil()
            logger.info("Fetch data from source %s in %s", urlTargetPfam, dirPath)
            fp = os.path.join(dirPath, fU.getFileName(urlTargetPfam))
            ok = fU.get(urlTargetPfam, fp)
            if not ok:
                fp = os.path.join(dirPath, fU.getFileName(urlTargetPfamFB))
                ok = fU.get(urlTargetPfamFB, fp)
                logger.info("Fetch data fallback fetch status is %r", ok)
            pfamD = self.__getPfamIndex(fp)
            ok = self.__mU.doExport(pfamDataPath, pfamD, fmt=fmt)
            logger.info("Caching %d in %s status %r", len(pfamD), pfamDataPath, ok)
            # ------
        #
        return pfamD

    def __getPfamIndex(self, filePath):
        """ Parse
        #
        """
        pfamD = {}
        encodingD = {"encoding": "ascii"} if sys.version_info[0] < 3 else {}
        rowL = self.__mU.doImport(filePath, fmt="tdd", rowFormat="list", **encodingD)
        for row in rowL:
            try:
                pfamId = row[0].strip().upper()
                idCode = row[3].strip()
                descr = row[4].strip()
                pfamD[pfamId] = descr + " (" + idCode + ")"
            except Exception:
                pass
        #
        return pfamD
