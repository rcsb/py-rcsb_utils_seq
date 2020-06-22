##
# File:    SiftsSummaryProvider.py
# Author:  J. Westbrook
# Date:    11-Dec-2018
#
# Updates:
#  06-Oct-2019 jdw add GO as part of abbreviated mappings
##
"""
Utilities to manage access to SIFTS summary mapping data.

"""
__docformat__ = "restructuredtext en"
__author__ = "John Westbrook"
__email__ = "jwest@rcsb.rutgers.edu"
__license__ = "Apache 2.0"

import logging
import os
import sys

from rcsb.utils.io.MarshalUtil import MarshalUtil
from rcsb.utils.seq.SeqAlign import SeqAlign, splitSeqAlignObjList

logger = logging.getLogger(__name__)


class SiftsSummaryProvider(object):
    """ Utilities to manage access to SIFTS summary mapping data.
    """

    def __init__(self, **kwargs):
        self.__ssD = self.__rebuildCache(**kwargs)

    def getEntryCount(self):
        return len(self.__ssD)

    def getEntries(self):
        return sorted(self.__ssD.keys())

    def getUniqueIdentifiers(self, idType="UNPID"):
        uL = []
        if idType not in ["UNPAL", "UNPID", "PFAMID", "GOID", "IPROID", "TAXID", "CATHID", "SCOPID", "ECID"]:
            return uL
        try:
            for _, aD in self.__ssD.items():
                for _, iD in aD.items():
                    if idType in iD:
                        uL.extend(iD[idType])
            return sorted(set(uL))
        except Exception:
            pass
        return uL

    def getEntryUniqueIdentifiers(self, entryIdList, idType="UNPID"):
        uL = []
        if idType not in ["UNPAL", "UNPID", "PFAMID", "GOID", "IPROID", "TAXID", "CATHID", "SCOPID", "ECID"]:
            return uL
        try:
            for entryId in entryIdList:
                if entryId in self.__ssD:
                    aD = self.__ssD[entryId]
                    for _, iD in aD.items():
                        if idType in iD:
                            uL.extend(iD[idType])
            return sorted(set(uL))
        except Exception:
            pass
        return uL

    def getAlignmentCount(self, entryId, authAsymId):
        num = 0
        try:
            num = len(self.__ssD[entryId][authAsymId]["UNPAL"])
        except Exception:
            pass
        return num

    def getAlignments(self, entryId, authAsymId):
        aL = []
        try:
            aL = self.__ssD[entryId][authAsymId]["UNPAL"]
        except Exception:
            pass
        return aL

    def getSeqAlignObjList(self, entryId, authAsymId):
        saoL = []
        try:
            # aL = self.__ssD[entryId][authAsymId]["UNPAL"]
            saoL = [SeqAlign("SIFTS", **sa) for sa in self.__ssD[entryId][authAsymId]["UNPAL"]]
        except Exception:
            pass
        return saoL

    def getLongestAlignments(self, entryId, authAsymIdL):
        """ Return the longest unique SIFTS alignments for the input entity instances.

        Args:
            entryId (str): entry identifier
            authAsymIdL (list): list of author entity instance (chain) identifiers

        Returns:
            (dict): {(dbName,dbAccession): [List of SeqAlign() objects with greatest coverage], }
        """
        retD = {}
        seqAlignObjL = []
        for authAsymId in authAsymIdL:
            seqAlignObjL.extend([SeqAlign("SIFTS", **sa) for sa in self.getIdentifiers(entryId, authAsymId, idType="UNPAL")])
        #
        alRefD = {}
        for seqAlignObj in seqAlignObjL:
            alRefD.setdefault((seqAlignObj.getDbName(), seqAlignObj.getDbAccession()), []).append(seqAlignObj)
        #
        # Get the longest overlapping entity region of each ref alignment -
        for (dbName, dbAcc), aL in alRefD.items():
            alGrpD = splitSeqAlignObjList(aL)
            logger.debug("SIFTS -> entryId %s dbName %r dbAcc %r alGrpD %r", entryId, dbName, dbAcc, alGrpD)
            for _, grpAlignL in alGrpD.items():
                lenL = [seqAlignObj.getEntityAlignLength() for seqAlignObj in grpAlignL]
                idxMax = lenL.index(max(lenL))
                retD.setdefault((dbName, dbAcc), []).append(grpAlignL[idxMax])
        return retD

    def getIdentifiers(self, entryId, authAsymId, idType=None):
        aL = []
        try:
            if idType in ["UNPAL", "UNPID", "PFAMID", "GOID", "IPROID", "TAXID", "CATHID", "SCOPID", "ECID"]:
                aL = self.__ssD[entryId][authAsymId][idType]
                logger.debug("sifts returns entryId %r authasymid %r idtype %r result %r", entryId, authAsymId, idType, aL)
            else:
                logger.error("Unsupported SIFTS idType %r", idType)
        except Exception as e:
            logger.debug("Failing with %s", str(e))
        return aL

    def getTaxIds(self, entryId, authAsymId):
        tL = []
        try:
            tL = self.__ssD[entryId][authAsymId]["TAXID"]
        except Exception:
            pass
        return tL

    def testCache(self):
        logger.info("SIFTS entry length %d", self.getEntryCount())
        if self.getEntryCount() > 140000:
            return True
        return False

    def __rebuildCache(self, **kwargs):
        mU = MarshalUtil()
        # source directory path
        srcDirPath = kwargs.get("srcDirPath", None)
        # cache details
        cacheKwargs = kwargs.get("cacheKwargs", {"fmt": "pickle"})
        useCache = kwargs.get("useCache", True)
        entrySaveLimit = kwargs.get("entrySaveLimit", None)
        abbreviated = str(kwargs.get("abbreviated", "TEST")).upper()
        #
        cacheDirPath = kwargs.get("cacheDirPath", None)
        pyVersion = sys.version_info[0]
        ext = "pic" if cacheKwargs["fmt"] == "pickle" else "json"
        saveFilePath = os.path.join(cacheDirPath, "sifts-summary-py%s.%s" % (str(pyVersion), ext))
        #
        ssD = {}
        try:
            if useCache and os.access(saveFilePath, os.R_OK):
                ssD = mU.doImport(saveFilePath, **cacheKwargs)
            else:
                if not srcDirPath:
                    logger.error("Missing SIFTS source path details")
                    return ssD
                ssD = self.__getSummaryMapping(srcDirPath, abbreviated=abbreviated)
                if entrySaveLimit:
                    ssD = {k: ssD[k] for k in list(ssD.keys())[:entrySaveLimit]}
                mU.mkdir(cacheDirPath)
                ok = mU.doExport(saveFilePath, ssD, **cacheKwargs)
                logger.debug("Saving SIFTS summary serialized data file %s (%d) status %r", saveFilePath, len(ssD), ok)
        except Exception as e:
            logger.exception("Failing with %s", str(e))
        return ssD

    def __getSummaryMapping(self, siftsSummaryDirPath, abbreviated="PROD"):
        """
        """

        uSeqD = self.__getUniprotChainMapping(siftsSummaryDirPath, "pdb_chain_uniprot.csv.gz")
        # _, uSeqD = self.__getUniprotChainMapping(siftsSummaryDirPath, "uniprot_segments_observed.csv.gz")
        logger.debug("uSeqD %d", len(uSeqD))
        for entryId, eD in uSeqD.items():
            for chainId, _ in eD.items():
                uSeqD[entryId][chainId]["UNPID"] = (
                    sorted(set(uSeqD[entryId][chainId]["UNPID"])) if ((entryId in uSeqD) and (chainId in uSeqD[entryId]) and ("UNPID" in uSeqD[entryId][chainId])) else []
                )

        #
        #
        tD = self.__getPfamChainMapping(siftsSummaryDirPath, "pdb_chain_pfam.csv.gz")
        logger.info("SIFTS PFAM mapping length %d", len(tD))
        for entryId, eD in uSeqD.items():
            for chainId, _ in eD.items():
                uSeqD[entryId][chainId]["PFAMID"] = sorted(set(tD[entryId][chainId])) if entryId in tD and chainId in tD[entryId] else []
        #
        if abbreviated == "TEST":
            return uSeqD
        #
        tD = self.__getInterProChainMapping(siftsSummaryDirPath, "pdb_chain_interpro.csv.gz")
        logger.info("SIFTS InterPro mapping length %d", len(tD))
        #
        for entryId, eD in uSeqD.items():
            for chainId, _ in eD.items():
                uSeqD[entryId][chainId]["IPROID"] = sorted(set(tD[entryId][chainId])) if entryId in tD and chainId in tD[entryId] else []

        #
        tD = self.__getGoIdChainMapping(siftsSummaryDirPath, "pdb_chain_go.csv.gz")
        logger.info("SIFTS GO mapping length %d", len(tD))
        #
        for entryId, eD in uSeqD.items():
            for chainId, _ in eD.items():
                uSeqD[entryId][chainId]["GOID"] = sorted(set(tD[entryId][chainId])) if entryId in tD and chainId in tD[entryId] else []
        #
        if abbreviated == "PROD":
            return uSeqD
        # --------------------
        tD = self.__getTaxonomnyChainMapping(siftsSummaryDirPath, "pdb_chain_taxonomy.csv.gz")
        logger.info("SIFTS Taxonomy mapping length %d", len(tD))
        for entryId, eD in uSeqD.items():
            for chainId, _ in eD.items():
                uSeqD[entryId][chainId]["TAXID"] = sorted(list(tD[entryId][chainId].keys())) if entryId in tD and chainId in tD[entryId] else []
        #
        tD = self.__getCathChainMapping(siftsSummaryDirPath, "pdb_chain_cath_uniprot.csv.gz")
        logger.info("SIFTS CATH mappinglength %d", len(tD))
        for entryId, eD in uSeqD.items():
            for chainId, _ in eD.items():
                uSeqD[entryId][chainId]["CATHID"] = tD[entryId][chainId] if entryId in tD and chainId in tD[entryId] else []

        tD = self.__getScopChainMapping(siftsSummaryDirPath, "pdb_chain_scop_uniprot.csv.gz")
        logger.info("SIFTS SCOP mapping length %d", len(tD))
        for entryId, eD in uSeqD.items():
            for chainId, _ in eD.items():
                uSeqD[entryId][chainId]["SCOPID"] = tD[entryId][chainId] if entryId in tD and chainId in tD[entryId] else []
        #

        tD = self.__getEnzymeChainMapping(siftsSummaryDirPath, "pdb_chain_enzyme.csv.gz")
        logger.info("SIFTS EC mapping length %d", len(tD))
        for entryId, eD in uSeqD.items():
            for chainId, _ in eD.items():
                uSeqD[entryId][chainId]["ECID"] = tD[entryId][chainId] if entryId in tD and chainId in tD[entryId] else []
        #
        return uSeqD

    def __getGeneOntologyChainMapping(self, siftsSummaryDirPath, csvFileName):
        """ Integrated SIFTS summary instance-level GO mapping data.

        PDB,CHAIN,SP_PRIMARY,WITH_STRING,EVIDENCE,GO_ID
        101m,A,IPRO,InterPro:IPR000971,IEA,GO:0020037
        101m,A,IPRO,InterPro:IPR002335,IEA,GO:0015671
        101m,A,IPRO,InterPro:IPR002335,IEA,GO:0019825
        101m,A,IPRO,InterPro:IPR002335,IEA,GO:0020037
        101m,A,IPRO,InterPro:IPR012292,IEA,GO:0019825

        """
        fp = os.path.join(siftsSummaryDirPath, csvFileName)
        rowDL = self.__readSiftsSummaryFile(fp)
        logger.info("Length of SIFTS summary file %s %d", csvFileName, len(rowDL))
        logger.debug("CSV keys: %r", list(rowDL[0].items()))
        tD = {}
        for rowD in rowDL:
            entryId = rowD["PDB"]
            chainId = rowD["CHAIN"]
            prov = rowD["SP_PRIMARY"]
            evidenceCode = rowD["EVIDENCE"]
            goId = rowD["GO_ID"]
            dD = {"GO_ID": goId, "EV": evidenceCode, "PROV": prov}
            tD.setdefault(entryId.upper(), {}).setdefault(chainId, []).append(dD)
        logger.info("GO data for %d entries", len(tD))
        return tD

    def __getGoIdChainMapping(self, siftsSummaryDirPath, csvFileName):
        """ Integrated SIFTS summary instance-level GO mapping data.

        PDB,CHAIN,SP_PRIMARY,WITH_STRING,EVIDENCE,GO_ID
        101m,A,IPRO,InterPro:IPR000971,IEA,GO:0020037
        101m,A,IPRO,InterPro:IPR002335,IEA,GO:0015671
        101m,A,IPRO,InterPro:IPR002335,IEA,GO:0019825
        101m,A,IPRO,InterPro:IPR002335,IEA,GO:0020037
        101m,A,IPRO,InterPro:IPR012292,IEA,GO:0019825

        """
        fp = os.path.join(siftsSummaryDirPath, csvFileName)
        # read this in list format as the number fields per record is variable.
        rowLL = self.__readSiftsSummaryFile(fp, rowFormat="list")
        logger.info("Length of SIFTS summary file %s %d", csvFileName, len(rowLL))
        # logger.debug("CSV keys: %r", list(rowDL[0].items()))
        tD = {}
        okCount = 0
        errCount = 0
        for rowL in rowLL:
            entryId = rowL[0]
            chainId = rowL[1]
            goId = rowL[-1]
            if not goId.startswith("GO:"):
                logger.warning("Skipping bad GO record (%d) %s %s %r", len(rowL), entryId, chainId, goId)
                errCount += 1
                continue
            tD.setdefault(entryId.upper(), {}).setdefault(chainId, []).append(goId)
            okCount += 1

        logger.info("GO data for %d entries", len(tD))
        logger.info("GO records %d format errors %d", okCount, errCount)
        return tD

    def __getInterProChainMapping(self, siftsSummaryDirPath, csvFileName):
        """ Integrated SIFTS summary instance-level InterPro mapping data.

        """
        fp = os.path.join(siftsSummaryDirPath, csvFileName)
        rowDL = self.__readSiftsSummaryFile(fp)
        logger.info("Length of SIFTS summary file %s %d", csvFileName, len(rowDL))
        logger.debug("CSV keys: %r", list(rowDL[0].items()))
        tD = {}
        for rowD in rowDL:
            entryId = rowD["PDB"]
            chainId = rowD["CHAIN"]
            interProId = rowD["INTERPRO_ID"]
            if not interProId.startswith("IPR"):
                logger.warning("Skipping bad InterPro ID %s %s %r %r", entryId, chainId, interProId, rowD)
                continue
            tD.setdefault(entryId.upper(), {}).setdefault(chainId, []).append(interProId)
        logger.info("InterPro data for %d entries", len(tD))
        return tD

    def __getPfamChainMapping(self, siftsSummaryDirPath, csvFileName):
        """ Integrated SIFTS summary instance-level PFAM mapping data.

        """
        fp = os.path.join(siftsSummaryDirPath, csvFileName)
        rowDL = self.__readSiftsSummaryFile(fp)
        logger.info("Length of SIFTS summary file %s %d", csvFileName, len(rowDL))
        logger.debug("CSV keys: %r", list(rowDL[0].items()))
        tD = {}
        for rowD in rowDL:
            entryId = rowD["PDB"]
            chainId = rowD["CHAIN"]
            pfamId = rowD["PFAM_ID"]
            if not pfamId.startswith("PF"):
                logger.warning("Skipping bad pfam ID %s %s %r %r", entryId, chainId, pfamId, rowD)
                continue
            tD.setdefault(entryId.upper(), {}).setdefault(chainId, []).append(pfamId)
        logger.info("PFAM data for %d entries", len(tD))
        return tD

    def __getEnzymeChainMapping(self, siftsSummaryDirPath, csvFileName):
        """ Integrated SIFTS summary instance-level EC mapping data.

        """
        fp = os.path.join(siftsSummaryDirPath, csvFileName)
        rowDL = self.__readSiftsSummaryFile(fp)
        logger.info("Length of SIFTS summary file %s %d", csvFileName, len(rowDL))
        logger.debug("CSV keys: %r", list(rowDL[0].items()))
        tD = {}
        for rowD in rowDL:
            entryId = rowD["PDB"]
            chainId = rowD["CHAIN"]
            ecId = rowD["EC_NUMBER"]
            tD.setdefault(entryId.upper(), {}).setdefault(chainId, []).append(ecId)
        logger.info("EC data for %d entries", len(tD))
        return tD

    def __getCathChainMapping(self, siftsSummaryDirPath, csvFileName):
        """ Integrated SIFTS summary instance-level CATH mapping data.

        """
        fp = os.path.join(siftsSummaryDirPath, csvFileName)
        rowDL = self.__readSiftsSummaryFile(fp)
        logger.info("Length of SIFTS summary file %s %d", csvFileName, len(rowDL))
        logger.debug("CSV keys: %r", list(rowDL[0].items()))
        tD = {}
        for rowD in rowDL:
            entryId = rowD["PDB"]
            chainId = rowD["CHAIN"]
            cathId = rowD["CATH_ID"]
            tD.setdefault(entryId.upper(), {}).setdefault(chainId, []).append(cathId)
        logger.info("CATH data for %d entries", len(tD))
        return tD

    def __getScopChainMapping(self, siftsSummaryDirPath, csvFileName):
        """ Integrated SIFTS summary instance-level SCOP mapping data.

        """
        fp = os.path.join(siftsSummaryDirPath, csvFileName)
        rowDL = self.__readSiftsSummaryFile(fp)
        logger.info("Length of SIFTS summary file %s %d", csvFileName, len(rowDL))
        logger.debug("%r", list(rowDL[0].items()))
        tD = {}
        for rowD in rowDL:
            entryId = rowD["PDB"]
            chainId = rowD["CHAIN"]
            scopId = rowD["SCOP_ID"]
            #
            tD.setdefault(entryId.upper(), {}).setdefault(chainId, []).append(scopId)
        #
        logger.info("SCOP data for %d entries", len(tD))
        return tD

    def __getTaxonomnyChainMapping(self, siftsSummaryDirPath, csvFileName):
        """ Integrated SIFTS summary instance-level taxonomy mapping data.

        # Taxonomy
        [('PDB', '101m'), ('CHAIN', 'A'), ('TAX_ID', '9755'), ('SCIENTIFIC_NAME', 'PHYCD')]

        """
        fp = os.path.join(siftsSummaryDirPath, csvFileName)
        rowDL = self.__readSiftsSummaryFile(fp)
        logger.info("Length of SIFTS summary file %s %d", csvFileName, len(rowDL))
        logger.debug("%r", list(rowDL[0].items()))
        tD = {}
        for rowD in rowDL:
            entryId = rowD["PDB"]
            chainId = rowD["CHAIN"]
            taxId = rowD["TAX_ID"]
            tD.setdefault(entryId.upper(), {}).setdefault(chainId, {}).update({taxId: True})
        #
        logger.info("Taxonomy for %d entries", len(tD))
        return tD

    def __getUniprotChainMapping(self, siftsSummaryDirPath, csvFileName):
        """  Integrated SIFTS summary instance-level uniprot mapping data.
        """
        #
        #
        fp = os.path.join(siftsSummaryDirPath, csvFileName)
        rowDL = self.__readSiftsSummaryFile(fp)
        logger.info("Length of SIFTS UniProt summary file %s %d", csvFileName, len(rowDL))
        logger.debug("%r", list(rowDL[0].items()))
        uD = {}
        # uIdD = {}
        for rowD in rowDL:
            entryId = rowD["PDB"]
            chainId = rowD["CHAIN"]
            unpId = rowD["SP_PRIMARY"]
            #
            entitySeqBeg = int(rowD["RES_BEG"]) if rowD["RES_BEG"].isdigit() else None
            entitySeqEnd = int(rowD["RES_END"]) if rowD["RES_END"].isdigit() else None
            entityLength = entitySeqEnd - entitySeqBeg + 1
            # authSeqBeg = int(rowD["PDB_BEG"]) if rowD["PDB_BEG"].isdigit() else None
            # authSeqEnd = int(rowD["PDB_END"]) if rowD["PDB_END"].isdigit() else None
            unpSeqBeg = int(rowD["SP_BEG"]) if rowD["SP_BEG"].isdigit() else None
            unpSeqEnd = int(rowD["SP_END"]) if rowD["SP_END"].isdigit() else None
            # dD = {"UP": unpId, "BG": entitySeqBeg, "ND": entitySeqEnd, "AUBG": authSeqBeg, "AUND": authSeqEnd, "UBG": unpSeqBeg, "UND": unpSeqEnd}
            # dD = {"UP": unpId, "BG": entitySeqBeg, "UBG": unpSeqBeg, "LEN": entityLength}
            dD = {"UP": unpId, "BG": entitySeqBeg, "LEN": entityLength, "UBG": unpSeqBeg, "UND": unpSeqEnd}
            uD.setdefault(entryId.upper(), {}).setdefault(chainId, {}).setdefault("UNPAL", []).append(dD)
            uD.setdefault(entryId.upper(), {}).setdefault(chainId, {}).setdefault("UNPID", []).append(unpId)
            #
        logger.info("UniProt mapping for %d entries", len(uD))
        # -----
        return uD

    def __readSiftsSummaryFile(self, filePath, rowFormat="dict"):
        """ Read input SIFTS summary file and return a list of dictionaries.
        """
        try:
            mU = MarshalUtil()
            cL = mU.doImport(filePath, fmt="csv", rowFormat=rowFormat)
            logger.debug("Container list %d", len(cL))
        except Exception as e:
            logger.exception("Failing with %s", str(e))
        return cL
