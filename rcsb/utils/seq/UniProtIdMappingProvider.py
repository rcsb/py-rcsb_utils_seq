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
from rcsb.utils.io.StashUtil import StashUtil

logger = logging.getLogger(__name__)


class UniProtIdMappingProvider(SingletonClass):
    """Manage index of UniProt identifier mappings."""

    def __init__(self, cachePath, **kwargs):
        #
        self.__cachePath = cachePath
        self.__dirName = "uniprot-id-mapping"
        self.__kwargs = kwargs
        self.__nameL = []
        self.__mapD = {}
        self.__cacheFieldD = {}
        self.__nameLegacyL = []
        self.__mapLegacyD = {}
        self.__cacheFieldLegacyD = {}
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

    def clearCache(self):
        try:
            self.__nameL = []
            self.__mapD = {}
            self.__cacheFieldD = {}
            self.__nameLegacyL = []
            self.__mapLegacyD = {}
            self.__cacheFieldLegacyD = {}
            dirPath = os.path.join(self.__cachePath, self.__dirName)
            fU = FileUtil()
            return fU.remove(dirPath)
        except Exception:
            pass
        return False

    def clearRawCache(self):
        try:
            rawDirPath = os.path.join(self.__cachePath, self.__dirName + "-raw")
            fU = FileUtil()
            return fU.remove(rawDirPath)
        except Exception:
            pass
        return False

    def reload(self, useCache=True, useLegacy=False, fmt="tdd", **kwargs):
        # full data set url -
        # urlTarget = kwargs.get("urlTarget", "ftp://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/idmapping/idmapping.dat.gz")
        ok = True
        urlTarget = kwargs.get("urlTarget", "ftp://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/idmapping/idmapping_selected.tab.gz")
        urlTargetLegacy = kwargs.get("urlTargetLegacy", "ftp://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/idmapping/idmapping_selected.tab.2015_03.gz")
        mapNameL = kwargs.get("mapNames", ["NCBI-taxon"])
        #
        dirPath = os.path.join(self.__cachePath, self.__dirName)
        rawDirPath = os.path.join(self.__cachePath, self.__dirName + "-raw")
        #
        if useLegacy and not (self.__mapLegacyD and len(self.__mapLegacyD) > 1000):
            logger.info("Reloading legacy mapping data (%r) from %r", useLegacy, urlTargetLegacy)
            self.__nameLegacyL, self.__mapLegacyD = self.__rebuildCache(urlTargetLegacy, mapNameL, dirPath, rawDirPath, fmt=fmt, useCache=useCache)
            self.__cacheFieldLegacyD = {k: ii for ii, k in enumerate(self.__nameLegacyL)}
            ok = self.__mapLegacyD and len(self.__mapLegacyD) > 1000
        elif not (self.__mapD and len(self.__mapD) > 1000):
            self.__nameL, self.__mapD = self.__rebuildCache(urlTarget, mapNameL, dirPath, rawDirPath, fmt=fmt, useCache=useCache)
            self.__cacheFieldD = {k: ii for ii, k in enumerate(self.__nameL)}
            ok = self.__mapD and len(self.__mapD) > 1000
        return ok

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
        logger.info("Length of UniProt ID mapping for %r %d", self.__nameL, len(self.__mapD))
        return (len(self.__mapD) > 1000) and (len(self.__nameL) >= 1)

    def __rebuildCache(self, targetUrl, mapNameL, outDirPath, rawDirPath, fmt="pickle", useCache=True):
        """Fetch the UniProt selected id mapping resource file and extract
        UniProt Acc to  'mapIndex' mapping. Serialize the mapping as required.

        Args:
            targetUrl (str): source URL of the remote index file
            mapNameL (list): list of key mapping names to extract from the index
            outDirPath (str): directory path for raw and processed mapping files
            fmt (str, optional): output format (pickle|json) . Defaults to "pickle".
            useCache (bool, optional): use cached files. Defaults to True.

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
                logger.info("Reading cached serialized file %r", idMapPath)
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
                idPath = os.path.join(rawDirPath, fileU.getFileName(targetUrl))
                if not fileU.exists(idPath):
                    logger.info("Fetching selected UniProt idmapping data from %r in %r", targetUrl, outDirPath)
                    ok = fileU.get(targetUrl, idPath)
                    if not ok:
                        logger.error("Failed to downlowd %r", targetUrl)
                        return oD
                else:
                    logger.info("Using cached mapping file %r", idPath)
                # ---
                ioU = IoUtil()
                if fmt in ["pickle", "json"]:
                    if len(mapNameL) == 1:
                        for row in ioU.deserializeCsvIter(idPath, delimiter="\t", rowFormat="list", encodingErrors="ignore"):
                            oD[row[0]] = str(row[self.__mapRecordD[mapNameL[0]] - 1])
                    else:
                        for row in ioU.deserializeCsvIter(idPath, delimiter="\t", rowFormat="list", encodingErrors="ignore"):
                            for mapName in mapNameL:
                                oD.setdefault(row[0], []).append(str(row[self.__mapRecordD[mapName] - 1]))
                    logger.info("Writing serialized mapping file %r", idMapPath)
                    ok = mU.doExport(idMapPath, {"idNameList": mapNameL, "uniprotMapD": oD}, fmt=fmt)
                elif fmt == "tdd":
                    #
                    logger.info("Writing serialized mapping file %r", idMapPath)
                    fU = FileUtil()
                    fU.mkdirForFile(idMapPath)
                    colNameL = []
                    colNameL.append("UniProtId")
                    colNameL.extend(mapNameL)
                    with open(idMapPath, "w") as ofh:
                        ofh.write("%s\n" % "\t".join(colNameL))
                        if len(mapNameL) == 1:
                            idx = self.__mapRecordD[mapNameL[0]] - 1
                            for row in ioU.deserializeCsvIter(idPath, delimiter="\t", rowFormat="list", encodingErrors="ignore"):
                                ofh.write("%s\t%s\n" % (row[0], row[idx]))
                        else:
                            idxL = [0]
                            idxL.extend([self.__mapRecordD[mapName] - 1 for mapName in mapNameL])
                            for row in ioU.deserializeCsvIter(idPath, delimiter="\t", rowFormat="list", encodingErrors="ignore"):
                                ofh.write("%s\n" % "\t".join([str(row[idx]) for idx in idxL]))
                            #
                    nL, oD = self.__rebuildCache(targetUrl, mapNameL, outDirPath, rawDirPath, fmt=fmt, useCache=True)
                    ok = True if nL and oD else False
            logger.info("Completed reload (%r) at %s (%.4f seconds)", ok, time.strftime("%Y %m %d %H:%M:%S", time.localtime()), time.time() - startTime)
        except Exception as e:
            logger.exception("Failing with %s", str(e))
        #
        return nL, oD

    def restore(self, cfgOb, configName):
        ok = False
        try:
            startTime = time.time()
            url = cfgOb.get("STASH_SERVER_URL", sectionName=configName)
            userName = cfgOb.get("_STASH_AUTH_USERNAME", sectionName=configName)
            password = cfgOb.get("_STASH_AUTH_PASSWORD", sectionName=configName)
            basePath = cfgOb.get("_STASH_SERVER_BASE_PATH", sectionName=configName)
            ok = self.__fromStash(url, basePath, userName=userName, password=password)
            logger.info("Recovered UniProt Mapping file from stash (%r)", ok)
            if not ok:
                urlFallBack = cfgOb.get("STASH_SERVER_FALLBACK_URL", sectionName=configName)
                ok = self.__fromStash(urlFallBack, basePath, userName=userName, password=password)
                logger.info("Recovered UniProt Mapping file from fallback stash (%r)", ok)
            #
            logger.info("Completed recovery (%r) at %s (%.4f seconds)", ok, time.strftime("%Y %m %d %H:%M:%S", time.localtime()), time.time() - startTime)
        except Exception as e:
            logger.exception("Failing with %s", str(e))
        #
        return ok

    def backup(self, cfgOb, configName):
        ok1 = ok2 = False
        try:
            startTime = time.time()
            userName = cfgOb.get("_STASH_AUTH_USERNAME", sectionName=configName)
            password = cfgOb.get("_STASH_AUTH_PASSWORD", sectionName=configName)
            basePath = cfgOb.get("_STASH_SERVER_BASE_PATH", sectionName=configName)
            url = cfgOb.get("STASH_SERVER_URL", sectionName=configName)
            urlFallBack = cfgOb.get("STASH_SERVER_FALLBACK_URL", sectionName=configName)
            ok1 = self.__toStash(url, basePath, userName=userName, password=password)
            ok2 = self.__toStash(urlFallBack, basePath, userName=userName, password=password)
            logger.info("Completed backup (%r/%r) at %s (%.4f seconds)", ok1, ok2, time.strftime("%Y %m %d %H:%M:%S", time.localtime()), time.time() - startTime)
        except Exception as e:
            logger.exception("Failing with %s", str(e))
        return ok1 & ok2

    def __toStash(self, url, stashRemoteDirPath, userName=None, password=None, remoteStashPrefix=None):
        """Copy tar and gzipped bundled cache data to remote server/location.

        Args:
            url (str): server URL (e.g. sftp://hostname.domain) None for local host
            stashRemoteDirPath (str): path to target directory on remote server
            userName (str, optional): server username. Defaults to None.
            password (str, optional): server password. Defaults to None.
            remoteStashPrefix (str, optional): channel prefix. Defaults to None.

        Returns:
            (bool): True for success or False otherwise
        """
        ok = False
        try:
            stU = StashUtil(os.path.join(self.__cachePath, "stash"), self.__dirName)
            ok = stU.makeBundle(self.__cachePath, [self.__dirName])
            if ok:
                ok = stU.storeBundle(url, stashRemoteDirPath, remoteStashPrefix=remoteStashPrefix, userName=userName, password=password)
        except Exception as e:
            logger.error("Failing with url %r stashDirPath %r: %s", url, stashRemoteDirPath, str(e))
        return ok

    def __fromStash(self, url, stashRemoteDirPath, userName=None, password=None, remoteStashPrefix=None):
        """Restore local cache from a tar and gzipped bundle to fetched from a remote server/location.

        Args:
            url (str): server URL (e.g. sftp://hostname.domain) None for local host
            stashRemoteDirPath (str): path to target directory on remote server
            userName (str, optional): server username. Defaults to None.
            password (str, optional): server password. Defaults to None.
            remoteStashPrefix (str, optional): channel prefix. Defaults to None.

        Returns:
            (bool): True for success or False otherwise
        """
        ok = False
        try:
            stU = StashUtil(os.path.join(self.__cachePath, "stash"), self.__dirName)
            ok = stU.fetchBundle(self.__cachePath, url, stashRemoteDirPath, remoteStashPrefix=remoteStashPrefix, userName=userName, password=password)
        except Exception as e:
            logger.error("Failing with url %r stashDirPath %r: %s", url, stashRemoteDirPath, str(e))
        return ok
