##
# File:    PfamProvider.py
# Date:    18-Feb-2020
#
# Updates:
#   26-May-2021 jdw Add methods to fetch and deliver Pfam-PDB mappings
#   20-Sep-2023 dwp Use HTTPS instead of FTP for Pfam data
#    3-Oct-2023 dwp Use new Pfam mapping file, pdbmap.gz, in place of pdb_pfamA_reg.txt.gz which is no longer updated/supported
#    3-Jun-2024 dwp Use new Pfam mapping file, pdb_pfam_mapping.tsv.gz (from SIFTS flatfiles),
#                   in place of pdbmap.gz which will eventually be deprecated
##

import logging
import os
import sys
from rcsb.utils.io.FileUtil import FileUtil
from rcsb.utils.io.MarshalUtil import MarshalUtil
from rcsb.utils.io.StashableBase import StashableBase

logger = logging.getLogger(__name__)


class PfamProvider(StashableBase):
    """Manage an index of Pfam identifier to description mappings."""

    def __init__(self, **kwargs):
        urlTargetPfam = kwargs.get("urlTargetPfam", "https://ftp.ebi.ac.uk/pub/databases/Pfam/current_release/Pfam-A.clans.tsv.gz")
        urlTargetPfamFB = "https://github.com/rcsb/py-rcsb_exdb_assets/raw/master/fall_back/Pfam/Pfam-A.clans.tsv.gz"
        self.__version = "34.0"
        dirName = "pfam"
        cachePath = kwargs.get("cachePath", ".")
        dirPath = os.path.join(cachePath, dirName)
        super(PfamProvider, self).__init__(cachePath, [dirName])
        useCache = kwargs.get("useCache", True)
        #
        self.__mU = MarshalUtil(workPath=dirPath)
        self.__pfamD = self.__rebuildCache(urlTargetPfam, urlTargetPfamFB, dirPath, useCache)

        urlTargetMapPfam = kwargs.get("urlTargetMapPfam", "https://ftp.ebi.ac.uk/pub/databases/msd/sifts/flatfiles/tsv/pdb_pfam_mapping.tsv.gz")
        urlTargetMapPfamFB = "https://github.com/rcsb/py-rcsb_exdb_assets/raw/master/fall_back/Pfam/pdb_pfam_mapping.tsv.gz"
        self.__pfamMapD = self.__rebuildMappingCache(urlTargetMapPfam, urlTargetMapPfamFB, dirPath, useCache)

    def getVersion(self):
        return self.__version

    def getDescription(self, pfamId):
        """Return the description for the input Pfam identifier

        Args:
            pfamId (str): Pfam identifier

        Returns:
            str: text description of the Pfam domain
        """
        descr = None
        try:
            descr = self.__pfamD[pfamId]
        except Exception:
            pass
        return descr

    def getMapping(self, pdbId):
        """Return the list of Pfam domain assignments for the input PDB identifer along with
        residue level mapping information

        Args:
            pdbId (str): PDB identifier

        Returns:
            list: [{'pfamId': , 'authAsymId":  , 'authSeqBeg': , 'authSeqEnd': 'insertBeg': , 'insertEnd': }, {}, ]
        """
        mapL = []
        try:
            mapL = self.__pfamMapD[pdbId.upper()]
        except Exception:
            pass
        return mapL

    def testCache(self):
        # Check length ...
        logger.info("Length pfamD %d pfamMapD %d", len(self.__pfamD), len(self.__pfamMapD))
        return (len(self.__pfamD) > 19000) and (len(self.__pfamMapD) > 150000)

    #
    def __rebuildCache(self, urlTargetPfam, urlTargetPfamFB, dirPath, useCache):
        pfamD = {}
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
        elif not useCache:
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
        """Parse annotation classifications
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

    def __rebuildMappingCache(self, urlTargetPfam, urlTargetPfamFB, dirPath, useCache):
        fmt = "json"
        ext = fmt if fmt == "json" else "pic"
        pfamDataPath = os.path.join(dirPath, "pfam-mapping-data.%s" % ext)
        #
        logger.debug("Using cache data path %s", dirPath)
        self.__mU.mkdir(dirPath)
        #
        if useCache and self.__mU.exists(pfamDataPath):
            pfamD = self.__mU.doImport(pfamDataPath, fmt=fmt)
            logger.debug("Pfam mapping data length %d", len(pfamD))
        else:
            # ------
            fU = FileUtil()
            logger.info("Fetch data from source %s in %s", urlTargetPfam, dirPath)
            fp = os.path.join(dirPath, fU.getFileName(urlTargetPfam))
            ok = fU.get(urlTargetPfam, fp)
            pfamD = {}
            if ok:
                pfamD = self.__getPfamMapping(fp)
            if not (ok and pfamD):
                fp = os.path.join(dirPath, fU.getFileName(urlTargetPfamFB))
                ok = fU.get(urlTargetPfamFB, fp)
                logger.info("Fetch data fallback fetch status is %r", ok)
                pfamD = self.__getPfamMapping(fp)
            ok = self.__mU.doExport(pfamDataPath, pfamD, fmt=fmt)
            logger.info("Caching %d in %s status %r", len(pfamD), pfamDataPath, ok)
            # ------
        #
        return pfamD

    def __getPfamMapping(self, filePath):
        """Parse mapping data"""
        pFamMapD = {}
        encodingD = {"encoding": "ascii"} if sys.version_info[0] < 3 else {}
        dataL = self.__mU.doImport(filePath, fmt="csv", csvDelimiter="\t", **encodingD)
        for rD in dataL:
            try:
                pdbId = rD["PDB"].strip().upper()
                pfamId = rD["PFAM_ACCESSION"].strip().upper()
                authAsymId = rD["CHAIN"].strip().upper()
                #
                authSeqBegRaw = rD["AUTH_PDBRES_START"].strip()
                authSeqBeg = int(authSeqBegRaw) if authSeqBegRaw and authSeqBegRaw != "None" else None
                #
                insertBegRaw = rD["AUTH_PDBRES_START_INS_CODE"].strip()
                insertBeg = insertBegRaw if insertBegRaw and insertBegRaw != "None" else None
                #
                authSeqEndRaw = rD["AUTH_PDBRES_END"].strip()
                authSeqEnd = int(authSeqEndRaw) if authSeqEndRaw and authSeqEndRaw != "None" else None
                #
                insertEndRaw = rD["AUTH_PDBRES_END_INS_CODE"].strip()
                insertEnd = insertEndRaw if insertEndRaw and insertEndRaw != "None" else None
                #
                pFamMapD.setdefault(pdbId, []).append(
                    {
                        "pfamId": pfamId,
                        "authAsymId": authAsymId,
                        "authSeqBeg": authSeqBeg,
                        "authSeqEnd": authSeqEnd,
                        "insertBeg": insertBeg,
                        "insertEnd": insertEnd,
                    }
                )
            except Exception as e:
                logger.exception("Failing with %r %s", rD, str(e))
        #
        logger.info("Pfam mapping data for (%d) entries", len(pFamMapD))
        return pFamMapD
