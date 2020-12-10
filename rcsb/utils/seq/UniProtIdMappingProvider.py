##
# File:    UniProtIdMappingProvider.py
# Date:    30-Nov-2020
#
##

import logging
import os
import time

from rcsb.utils.io.FileUtil import FileUtil
from rcsb.utils.io.IoUtil import IoUtil
from rcsb.utils.io.MarshalUtil import MarshalUtil
from rcsb.utils.io.SingletonClass import SingletonClass

logger = logging.getLogger(__name__)


class UniProtIdMappingProvider(SingletonClass):
    """Manage index of UniProt identifier mappings."""

    def __init__(self, **kwargs):
        # urlTarget = kwargs.get("urlTarget", "ftp://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/idmapping/idmapping.dat.gz")
        urlTarget = kwargs.get("urlTarget", "ftp://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/idmapping/idmapping_selected.tab.gz")
        urlTargetLegacy = kwargs.get("urlTargetLegacy", "ftp://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/idmapping/idmapping_selected.tab.2015_03.gz")
        cachePath = kwargs.get("cachePath", ".")
        dirPath = os.path.join(cachePath, "uniprot-id-mapping")
        useCache = kwargs.get("useCache", True)
        mapNameL = kwargs.get("mapNames", ["NCBI-taxon"])
        fmt = kwargs.get("fmt", "pickle")
        maxLimit = kwargs.get("maxLimit", None)
        useLegacy = kwargs.get("useLegacy", False)
        #
        self.__mapRecordD = {
            "UniProtKB-AC": 1,
            "UniProtKB-ID": 2,
            "GeneID": 3,
            "RefSeq": 4,
            "GI": 5,
            "PDB": 6,
            "GO": 7,
            "UniRef100": 8,
            "UniRef90": 9,
            "UniRef50": 10,
            "UniParc": 11,
            "PIR": 12,
            "NCBI-taxon": 13,
            "MIM": 14,
            "UniGene": 15,
            "PubMed": 16,
            "EMBL": 17,
            "EMBL-CDS": 18,
            "Ensembl": 19,
            "Ensembl_TRS": 20,
            "Ensembl_PRO": 21,
            "Additional PubMed": 22,
        }
        self.__mU = MarshalUtil(workPath=dirPath)
        self.__nameL, self.__mapD = self.__rebuildCache(urlTarget, mapNameL, dirPath, fmt=fmt, useCache=useCache, maxLimit=maxLimit)
        self.__cacheFieldD = {k: ii for ii, k in enumerate(self.__nameL)}
        #
        if useLegacy:
            self.__nameLegacyL, self.__mapLegacyD = self.__rebuildCache(urlTargetLegacy, mapNameL, dirPath, fmt=fmt, useCache=useCache, maxLimit=maxLimit)
            self.__cacheFieldLegacyD = {k: ii for ii, k in enumerate(self.__nameLegacyL)}
        else:
            self.__nameLegacyL = []
            self.__mapLegacyD = {}
            self.__cacheFieldLegacyD = {}

    def getMappedId(self, unpId, mapName="NCBI-taxon"):
        mId = None
        try:
            if len(self.__nameL) == 1 and self.__nameL[0] == mapName:
                mId = self.__mapD[unpId]
            else:
                iRec = self.__cacheFieldD[mapName]
                mId = self.__mapD[unpId][iRec]
        except Exception:
            pass
        return mId

    def getMappedIdLegacy(self, unpId, mapName="NCBI-taxon"):
        mId = None
        try:
            if len(self.__nameL) == 1 and self.__nameL[0] == mapName:
                mId = self.__mapLegacyD[unpId]
            else:
                iRec = self.__cacheFieldLegacyD[mapName]
                mId = self.__mapLegacyD[unpId][iRec]
        except Exception:
            pass
        return mId

    def testCache(self):
        logger.info("Length UniProt mapping for %r %d", self.__nameL, len(self.__mapD))
        return len(self.__mapD) > 1000

    def __rebuildCache(self, targetUrl, mapNameL, outDirPath, fmt="pickle", useCache=True, maxLimit=None):
        """Fetch the UniProt selected id mapping resource file and extract
        UniProt Acc to  'mapIndex' mapping. Serialize the mapping as required.

        Args:
            mapIndex (int): index in the tab delimited mapping file (below)
            outDirPath (str): directory path for raw and processed mapping files
            mapFileName (str): file name for the target mapping file
            fmt (str, optional): output format (pickle|json) . Defaults to "pickle".
            useCache (bool, optional): use cached files. Defaults to True.
            maxLimit (int, optional): maximum number of records to process (default=None)

        Returns:
            dict: od[uniprotId] = mapped value

                idmapping_selected.tab

                1. UniProtKB-AC
                2. UniProtKB-ID
                3. GeneID (EntrezGene)
                4. RefSeq
                5. GI
                6. PDB
                7. GO
                8. UniRef100
                9. UniRef90
                10. UniRef50
                11. UniParc
                12. PIR
                13. NCBI-taxon
                14. MIM
                15. UniGene
                16. PubMed
                17. EMBL
                18. EMBL-CDS
                19. Ensembl
                20. Ensembl_TRS
                21. Ensembl_PRO
                22. Additional PubMed

        """
        startTime = time.time()
        nL = mapNameL
        oD = {}
        try:
            fileU = FileUtil()
            fExt = "pic" if fmt == "pickle" else "json"
            fExt = "tdd" if fmt == "tdd" else fExt
            fN, _ = os.path.splitext(fileU.getFileName(targetUrl))
            mapFileName = fN + "-map." + fExt
            idMapPath = os.path.join(outDirPath, mapFileName)
            mU = MarshalUtil()
            if useCache and mU.exists(idMapPath):
                logger.info("Using cached file %r", idMapPath)
                if fmt in ["pickle", "json"]:
                    tD = mU.doImport(idMapPath, fmt=fmt)
                    nL = list(set(tD["idNameList"]))
                    oD = tD["uniprotMapD"]
                    logger.info("keys %r", list(oD.keys())[:10])
                    logger.info("nL %r", nL)
                    ok = True
                elif fmt == "tdd":
                    ioU = IoUtil()
                    it = ioU.deserializeCsvIter(idMapPath, delimiter="\t", rowFormat="list", encodingErrors="ignore")
                    tL = next(it, [])
                    nL = tL[1:]
                    if len(nL) == 1:
                        for row in it:
                            oD[row[0]] = row[1]
                    else:
                        for row in it:
                            oD[row[0]] = row[1:]
                    ok = True
            else:
                idPath = os.path.join(outDirPath, fileU.getFileName(targetUrl))
                if not fileU.exists(idPath):
                    logger.info("Fetch selected idmapping data from %r in %r", targetUrl, outDirPath)
                    ok = fileU.get(targetUrl, idPath)
                    if not ok:
                        logger.error("Failed to downlowd %r", targetUrl)
                        return oD
                else:
                    logger.info("Using cached mapping file %r", idPath)
                # ---
                ioU = IoUtil()
                if fmt in ["pickle", "json"]:
                    iCount = 0
                    if len(mapNameL) == 1:
                        for row in ioU.deserializeCsvIter(idPath, delimiter="\t", rowFormat="list", encodingErrors="ignore"):
                            oD[row[0]] = str(row[self.__mapRecordD[mapNameL[0]] - 1])
                            iCount += 1
                            if maxLimit and iCount > maxLimit:
                                break
                    else:
                        for row in ioU.deserializeCsvIter(idPath, delimiter="\t", rowFormat="list", encodingErrors="ignore"):
                            for mapName in mapNameL:
                                oD.setdefault(row[0], []).append(str(row[self.__mapRecordD[mapName] - 1]))
                            iCount += 1
                            if maxLimit and iCount > maxLimit:
                                break
                    ok = mU.doExport(idMapPath, {"idNameList": mapNameL, "uniprotMapD": oD}, fmt=fmt)
                elif fmt == "tdd":
                    #
                    colNameL = []
                    colNameL.append("UniProtId")
                    colNameL.extend(mapNameL)
                    oL = []
                    oL.append("\t".join(colNameL))

                    iCount = 0
                    if len(mapNameL) == 1:
                        for row in ioU.deserializeCsvIter(idPath, delimiter="\t", rowFormat="list", encodingErrors="ignore"):
                            oL.append(row[0] + "\t" + row[self.__mapRecordD[mapNameL[0]] - 1])
                            iCount += 1
                            if maxLimit and iCount > maxLimit:
                                break
                    else:
                        for row in ioU.deserializeCsvIter(idPath, delimiter="\t", rowFormat="list", encodingErrors="ignore"):
                            tL = [row[0]]
                            for mapName in mapNameL:
                                tL.append(str(row[self.__mapRecordD[mapName] - 1]))
                            oL.append("\t".join(tL))
                            iCount += 1
                            if maxLimit and iCount > maxLimit:
                                break
                        #
                    ok = mU.doExport(idMapPath, oL, fmt="list")
                    oL = []
                    nL, oD = self.__rebuildCache(targetUrl, mapNameL, outDirPath, fmt=fmt, useCache=True, maxLimit=maxLimit)
            logger.info("Completed reload (%r) at %s (%.4f seconds)", ok, time.strftime("%Y %m %d %H:%M:%S", time.localtime()), time.time() - startTime)
        except Exception as e:
            logger.exception("Failing with %s", str(e))
        #
        return nL, oD
