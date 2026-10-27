##
#  File:           PlasticEnzymeAnnotationProvider.py
#  Date:           7-Oct-2026
#
#  Updates:
#  20-Oct-2026 dwp Use the PlasticDB FASTA file (instead of the TSV) as the PlasticDB data source; parse all PlasticDB
#                  FASTA header fields; read FASTA files with MarshalUtil; skip writing source data with too few
#                  records; add reload/testCache, buildAnnotationList (list of annotations per entity) and accessors;
#                  sort substrate lists; don't overwrite the annotation file when the search results are empty or
#                  yield too few annotated entities; send a 'RCSB-PDB-ETL' User-Agent with source data requests;
#                  log failed source fetches (e.g., HTTP 403) as errors with the source URL
#  26-Oct-2026 dwp Ignore (and log) cache files that don't load as a dictionary
#  27-Oct-2026 dwp Add a request timeout for source data fetches (default 30 seconds), so a stalled connection
#                  fails instead of hanging; download source files to a temporary file and keep them (pazy_proteins.json,
#                  PlasticDB.fasta) only once they parse with enough records, so the raw files are part of the stashed
#                  backup; if either source can't be fetched, restore the stashed backup and write the freshly fetched
#                  source over it; restore() uses a valid local cache (testCache()) instead of the stash; backup() first
#                  adds any missing raw source files from the current stashed backup, and skips the backup without them
##
"""
Accessors and generators for plastic-degrading enzyme annotation data (sourced from PAZy and PlasticDB).

Workflow:
    1. fetchPlasticEnzymeData() - fetch PAZy (JSON) and PlasticDB (FASTA) data, keep the raw files (pazy_proteins.json
                                  and PlasticDB.fasta) and standardize each (pazy-select-data.json and
                                  plasticdb-select-data.json), with fallback to the stashed backup
    2. exportFasta()            - export the combined enzyme sequences (plastic-enzyme-sequences.fa)
    3. (external)               - sequence search of the enzyme sequences against PDB polymer entities
    4. buildAnnotationList()    - build the polymer entity annotations from the search results
                                  (plastic-enzyme-annotation-data.json)
    5. getAnnotation()          - annotation lookup by polymer entity ID
"""

import datetime
import logging
import os.path
import re
import shutil
import tempfile
import time

from rcsb.utils.io.FileUtil import FileUtil
from rcsb.utils.io.MarshalUtil import MarshalUtil
from rcsb.utils.io.StashableBase import StashableBase

logger = logging.getLogger(__name__)


class PlasticEnzymeAnnotationProvider(StashableBase):
    """Accessors and generators for plastic-degrading enzyme annotation data."""

    def __init__(self, **kwargs):
        #
        self.__cachePath = kwargs.get("cachePath", ".")
        self.__dirName = "plastic-enzyme-annotations"
        super(PlasticEnzymeAnnotationProvider, self).__init__(self.__cachePath, [self.__dirName])
        self.__dirPath = os.path.join(self.__cachePath, self.__dirName)
        # Raw source data files, kept in the cache directory (and so in the stashed backup)
        self.__rawFileNameD = {"PAZy": "pazy_proteins.json", "PlasticDB": "PlasticDB.fasta"}
        #
        # Source locators (remote URLs by default; local file paths are also accepted by FileUtil.get())
        self.__pazyUrl = kwargs.get("pazyUrl", "https://api.pazy.eu/api/proteins/json/?format=json")
        self.__plasticDBUrl = kwargs.get("plasticDBUrl", "https://plasticdb.org/static/PlasticDB.fasta")
        userAgent = kwargs.get("userAgent", "RCSB-PDB-ETL")
        self.__headers = {"User-Agent": userAgent} if userAgent else None
        # Seconds without a response (connect, TLS handshake or between bytes) before a source fetch attempt fails
        self.__timeout = kwargs.get("timeout", 30)
        #
        self.__mU = MarshalUtil(workPath=self.__dirPath)
        self.__fU = FileUtil(timeout=self.__timeout)
        self.__pazyD, self.__plasticDBD, self.__annotD = self.__reload(**kwargs)

    def testCache(self, minCount=100):
        """Check that the PAZy and PlasticDB source data and the annotation data are each loaded with more than minCount records."""
        ok = self.testSourceCache(minCount=minCount)
        annotCount = len(self.__annotD.get("annotations", {})) if self.__annotD else 0
        logger.info("Plastic enzyme annotation count %d (minCount %d)", annotCount, minCount)
        return ok and annotCount > minCount

    def testSourceCache(self, minCount=100):
        """Check that the PAZy and PlasticDB source data are each loaded with more than minCount records (required for
        exportFasta() and buildAnnotationList())."""
        pazyCount = len(self.__pazyD.get("data", {})) if self.__pazyD else 0
        plasticDBCount = len(self.__plasticDBD.get("data", {})) if self.__plasticDBD else 0
        logger.info("Plastic enzyme source data counts PAZy %d PlasticDB %d (minCount %d)", pazyCount, plasticDBCount, minCount)
        return pazyCount > minCount and plasticDBCount > minCount

    def reload(self):
        """Reload the source and annotation data from the cache path.

        Returns:
            bool: testCache() status
        """
        self.__pazyD, self.__plasticDBD, self.__annotD = self.__reload(useCache=True)
        return self.testCache()

    def getVersion(self):
        return self.__annotD.get("version", None) if self.__annotD else None

    def hasAnnotation(self, rcsbEntityId):
        return bool(self.__annotD) and rcsbEntityId in self.__annotD.get("annotations", {})

    def getAnnotation(self, rcsbEntityId):
        """Get the plastic enzyme annotations for the input polymer entity ID (e.g., "7DZT_1" or "pdb_00007dzt_1").

        Returns:
            list: annotation dictionaries (best match first), or None if the entity has no annotations
        """
        try:
            return self.__annotD["annotations"][rcsbEntityId]
        except Exception:
            pass
        return None

    def getEnzyme(self, sourceDb, enzymeId):
        """Get the standardized source record for the input source database ("PAZy" or "PlasticDB") and enzyme ID."""
        try:
            dD = self.__pazyD if sourceDb == "PAZy" else self.__plasticDBD if sourceDb == "PlasticDB" else {}
            return dD["data"][enzymeId]
        except Exception:
            pass
        return None

    def __getPazyDataPath(self):
        return os.path.join(self.__dirPath, "pazy-select-data.json")

    def __getPlasticDBDataPath(self):
        return os.path.join(self.__dirPath, "plasticdb-select-data.json")

    def __getSelectDataPath(self, sourceDb):
        return self.__getPazyDataPath() if sourceDb == "PAZy" else self.__getPlasticDBDataPath()

    def __getRawPath(self, sourceDb):
        return os.path.join(self.__dirPath, self.__rawFileNameD[sourceDb])

    def __getAnnotationDataPath(self):
        return os.path.join(self.__dirPath, "plastic-enzyme-annotation-data.json")

    def __getFastaPath(self):
        return os.path.join(self.__dirPath, "plastic-enzyme-sequences.fa")

    def __reload(self, **kwargs):
        startTime = time.time()
        useCache = kwargs.get("useCache", True)
        self.__fU.mkdir(self.__dirPath)
        dataL = [{}, {}, {}]
        if useCache:
            for ii, dataPath in enumerate([self.__getPazyDataPath(), self.__getPlasticDBDataPath(), self.__getAnnotationDataPath()]):
                if self.__mU.exists(dataPath):
                    dataD = self.__mU.doImport(dataPath, fmt="json")
                    if isinstance(dataD, dict):
                        dataL[ii] = dataD
                    else:
                        logger.warning("Unable to load plastic enzyme cache file %r", dataPath)
                else:
                    logger.info("Plastic enzyme cache file not present %r", dataPath)
        logger.info("Completed reload (useCache %r) at %s (%.4f seconds)", useCache, time.strftime("%Y %m %d %H:%M:%S", time.localtime()), time.time() - startTime)
        return dataL[0], dataL[1], dataL[2]

    def fetchPlasticEnzymeData(self, minCount=100, cfgOb=None, configName=None, remotePrefix=None):
        """Fetch the PAZy and PlasticDB raw source data and restructure each into a standardized JSON file
        (pazy-select-data.json and plasticdb-select-data.json) in the cache path. The raw source files
        (pazy_proteins.json and PlasticDB.fasta) are kept there too, so they're part of the stashed backup.

        If either source can't be fetched, the stashed backup is restored (when cfgOb and configName are given), and
        the freshly fetched raw and standardized data for the other source are then written over the restored copies,
        so the newest data is kept. Use testSourceCache() to check that source data is available afterwards.

        Args:
            minCount (int, optional): minimum number of records required to accept the fetched data for a source. Defaults to 100.
            cfgOb (obj, optional): configuration object (ConfigUtil()) for restoring the stashed backup. Defaults to None (no restore).
            configName (str, optional): configuration section name for restoring the stashed backup. Defaults to None.
            remotePrefix (str, optional): stash channel prefix. Defaults to None.

        Returns:
            bool: True if both sources were freshly fetched or False otherwise
        """
        startTime = time.time()
        freshL = []
        failedL = []
        for sourceDb, fetchFunc in [("PAZy", self.__fetchPazyData), ("PlasticDB", self.__fetchPlasticDBData)]:
            dataD, downloadPath = fetchFunc(minCount=minCount)
            if dataD is None:
                failedL.append(sourceDb)
            else:
                freshL.append((sourceDb, dataD, downloadPath))
        #
        if failedL:
            if cfgOb and configName:
                logger.warning("Fetching %s data failed - restoring the stashed backup", " and ".join(failedL))
                okR = super(PlasticEnzymeAnnotationProvider, self).restore(cfgOb, configName, remotePrefix=remotePrefix, useStash=True, useGit=False)
                logger.info("Plastic enzyme stash restore status (%r)", okR)
            else:
                logger.error("Fetching %s data failed, and no stash configuration given for a restore", " and ".join(failedL))
        #
        # Keep the freshly fetched raw files and write their standardized data (after any stash restore, so they replace the restored copies)
        for sourceDb, dataD, downloadPath in freshL:
            os.replace(downloadPath, self.__getRawPath(sourceDb))
            self.__exportSelectData(self.__getSelectDataPath(sourceDb), dataD, sourceDb)
        self.__pazyD, self.__plasticDBD, _ = self.__reload(useCache=True)
        #
        ok = not failedL
        logger.info(
            "Completed plastic enzyme source data fetch (%r) fetched %r failed %r at %s (%.4f seconds)",
            ok,
            [tup[0] for tup in freshL],
            failedL,
            time.strftime("%Y %m %d %H:%M:%S", time.localtime()),
            time.time() - startTime,
        )
        return ok

    def restore(self, cfgOb, configName, remotePrefix=None, useStash=True, useGit=False):
        """Restore the plastic enzyme cache directory from the stashed backup, unless the local cache already holds
        valid data (testCache()), which is then used as is.

        Returns:
            bool: True if the local cache is valid or the restore succeeds, or False otherwise
        """
        self.__pazyD, self.__plasticDBD, self.__annotD = self.__reload(useCache=True)
        if self.testCache():
            logger.info("Using the valid local plastic enzyme cache in %r (no stash restore)", self.__dirPath)
            return True
        logger.info("Local plastic enzyme cache in %r is missing or incomplete - restoring the stashed backup", self.__dirPath)
        ok = super(PlasticEnzymeAnnotationProvider, self).restore(cfgOb, configName, remotePrefix=remotePrefix, useStash=useStash, useGit=useGit)
        self.__pazyD, self.__plasticDBD, self.__annotD = self.__reload(useCache=True)
        return ok

    def backup(self, cfgOb, configName, remotePrefix=None, useStash=True, useGit=False, backupToFallback=False):
        """Back up the plastic enzyme cache directory to the stash, first adding any raw source file that's missing
        locally (pazy_proteins.json or PlasticDB.fasta) from the current stashed backup, so that the new backup always
        includes both raw files.

        Returns:
            bool: True for success, or False otherwise (including when a raw source file can't be found, and the backup is skipped)
        """
        missingL = [sourceDb for sourceDb in self.__rawFileNameD if not self.__mU.exists(self.__getRawPath(sourceDb))]
        if missingL:
            logger.warning("Raw %s source files missing in %r - adding them from the stashed backup", " and ".join(missingL), self.__dirPath)
            tmpPath = tempfile.mkdtemp(prefix="plastic-enzyme-restore-", dir=self.__cachePath)
            try:
                # Restore the stashed backup into a separate directory, so the local data isn't overwritten
                tP = PlasticEnzymeAnnotationProvider(cachePath=tmpPath, useCache=False)
                StashableBase.restore(tP, cfgOb, configName, remotePrefix=remotePrefix, useStash=useStash, useGit=useGit)
                for sourceDb in list(missingL):
                    srcPath = os.path.join(tmpPath, self.__dirName, self.__rawFileNameD[sourceDb])
                    if self.__mU.exists(srcPath):
                        shutil.copyfile(srcPath, self.__getRawPath(sourceDb))
                        missingL.remove(sourceDb)
            except Exception as e:
                logger.exception("Failing to add raw source files from the stashed backup with %s", str(e))
            finally:
                shutil.rmtree(tmpPath, ignore_errors=True)
        if missingL:
            logger.error("Raw %s source files are missing locally and in the stashed backup - skipping the backup", " and ".join(missingL))
            return False
        return super(PlasticEnzymeAnnotationProvider, self).backup(cfgOb, configName, remotePrefix=remotePrefix, useStash=useStash, useGit=useGit, backupToFallback=backupToFallback)

    def __getSourceFile(self, sourceDb, url):
        """Download a raw source data file to a temporary file in the cache directory (removed on failure).

        Returns:
            str: path to the downloaded file, or None if it couldn't be fetched
        """
        downloadPath = self.__getRawPath(sourceDb) + ".download"
        logger.info("Fetching %s data from %r to %r", sourceDb, url, downloadPath)
        if not self.__fU.get(url, downloadPath, headers=self.__headers):
            self.__fU.remove(downloadPath)
            logger.error("Fetching %s data from %r failed (site down or blocking requests, e.g., HTTP 403 - see the HTTP status above)", sourceDb, url)
            return None
        return downloadPath

    def exportFasta(self, fastaPath=None):
        """Export a combined (PAZy and PlasticDB) plastic enzyme sequence fasta file (plastic-enzyme-sequences.fa)
        from the standardized source data.

        Each fasta header has the form: >{enzymeId}|enzymeId|{sourceDb}|sourceDb
        (e.g., >111|enzymeId|PAZy|sourceDb or >00003|enzymeId|PlasticDB|sourceDb).

        Args:
            fastaPath (str, optional): fasta output file path. Defaults to plastic-enzyme-sequences.fa in the cache path

        Returns:
            bool: True for success or False otherwise
        """
        fastaPath = fastaPath if fastaPath else self.__getFastaPath()
        ok = False
        try:
            if not (self.__pazyD and self.__plasticDBD):
                logger.error("Missing plastic enzyme source data (run fetchPlasticEnzymeData first)")
                return False
            sD = {}
            for qD in [self.__pazyD, self.__plasticDBD]:
                for enzymeId, tD in qD["data"].items():
                    if not tD.get("sequence"):
                        continue
                    # Dictionary order determines the fasta header: {value}|{key}|{value}|{key}
                    cD = {"sequence": tD["sequence"], "enzymeId": enzymeId, "sourceDb": tD["source_db"]}
                    cId = "|".join([str(enzymeId), "enzymeId", str(tD["source_db"]), "sourceDb"])
                    sD[cId] = cD
            #
            ok = self.__mU.doExport(fastaPath, sD, fmt="fasta", makeComment=True)
            logger.info("Export plastic enzyme fasta (%d) status %r", len(sD), ok)
        except Exception as e:
            logger.exception("Failing to export plastic enzyme fasta with %s", str(e))
        return ok

    def buildAnnotationList(self, sequenceMatchFilePath, minCount=100):
        """Build the polymer entity annotations for the matching entities in the input sequence match file
        -> plastic-enzyme-annotation-data.json

        A polymer entity can match several enzymes (e.g., identical sequences in PAZy and PlasticDB), so each entity
        is assigned a list of annotations (one per matching enzyme), sorted by best match first.

        Args:
            sequenceMatchFilePath (str): sequence match output file path (e.g., plastic-enzyme-vs-pdbprent-filtered-results.json)
            minCount (int, optional): minimum number of PAZy and PlasticDB source records required, and minimum number of
                                      annotated entities required to replace the cached annotation file. Defaults to 100.

        Returns:
            bool: True for success or False otherwise
        """
        if not self.testSourceCache(minCount=minCount):
            logger.warning("Skipping build of plastic enzyme annotations because the source data is missing")
            return False
        mD = self.__mU.doImport(sequenceMatchFilePath, fmt="json") if self.__mU.exists(sequenceMatchFilePath) else None
        if not mD:
            logger.error("Missing or empty sequence match file %r - keeping the existing annotation data", sequenceMatchFilePath)
            return False
        #
        # {entityId: {(sourceDb, enzymeId): annotation dict}}
        eD = {}
        for queryId, matchDL in mD.items():
            qCmtD = self.__decodeComment(queryId)
            enzymeId = qCmtD.get("enzymeId")
            sourceDb = qCmtD.get("sourceDb")
            tD = self.getEnzyme(sourceDb, enzymeId)
            if not tD:
                logger.warning("Skipping matches for unknown enzyme %r", queryId)
                continue
            for matchD in matchDL:
                try:
                    tCmtD = self.__decodeComment(matchD["target"])
                    rcsbEntityId = tCmtD["entityId"]
                    # Split on the last underscore, so that both short (7DZT_1) and extended (pdb_00007dzt_1) IDs work
                    entryId, entityId = rcsbEntityId.rsplit("_", 1)
                    seqIdPct = float(matchD["sequenceIdentity"])
                    rD = {
                        "entry_id": entryId,
                        "entity_id": entityId,
                        "annotation_id": enzymeId,
                        "source_db": sourceDb,
                        "name": tD.get("name"),
                        "description": "Plastic degrading enzyme (sequence matched against %s)." % sourceDb,
                        "substrates": tD.get("substrates", []),
                        "seq_id_pct": seqIdPct,
                        "bit_score": float(matchD["bitScore"]),
                        "reference_scheme": "PDB entity",
                        "perfect_match": "Y" if seqIdPct == 100.0 else "N",
                    }
                    # Keep the best match for each (entity, enzyme) pair (e.g., an entity with multiple source parts)
                    ky = (sourceDb, enzymeId)
                    oD = eD.setdefault(rcsbEntityId, {})
                    if ky not in oD or (rD["seq_id_pct"], rD["bit_score"]) > (oD[ky]["seq_id_pct"], oD[ky]["bit_score"]):
                        oD[ky] = rD
                except Exception as e:
                    logger.warning("Unable to process match data for %r with exception %r", queryId, e)
        #
        aD = {}
        for rcsbEntityId, oD in eD.items():
            aD[rcsbEntityId] = sorted(oD.values(), key=lambda k: (-k["seq_id_pct"], -k["bit_score"], k["source_db"], k["annotation_id"]))
        #
        if len(aD) < minCount:
            logger.error("Too few annotated entities (%d < %d) - keeping the existing annotation data", len(aD), minCount)
            return False
        tS = datetime.datetime.now().isoformat()
        self.__annotD = {"version": "1.0", "created": tS, "annotations": aD}
        ok = self.__mU.doExport(self.__getAnnotationDataPath(), self.__annotD, fmt="json", indent=3)
        logger.info("Built plastic enzyme annotations (status %r) for %d entities (%d annotations)", ok, len(aD), sum(len(v) for v in aD.values()))
        return ok

    def __decodeComment(self, comment, separator="|"):
        dD = {}
        try:
            ti = iter(comment.split(separator))
            dD = {tup[1]: tup[0] for tup in zip(ti, ti)}
        except Exception:
            pass
        return dD

    def __fetchPazyData(self, minCount=100):
        """Fetch PAZy JSON data and restructure it (written to pazy-select-data.json by fetchPlasticEnzymeData())

        Returns:
            (dict, str): standardized records keyed by enzyme ID and the path to the downloaded raw file, or (None, None) if
                         the data couldn't be fetched or has too few records (the download is then removed)
        """
        downloadPath = self.__getSourceFile("PAZy", self.__pazyUrl)
        if not downloadPath:
            return None, None
        rawD = self.__mU.doImport(downloadPath, fmt="json")
        if not rawD:
            self.__fU.remove(downloadPath)
            logger.error("Unable to read PAZy data from %r", self.__pazyUrl)
            return None, None
        # The bulk download is a dictionary with "count" and "results" keys (but also accept a plain list)
        if isinstance(rawD, dict):
            entryL = rawD.get("results", [])
            if "count" in rawD and rawD["count"] != len(entryL):
                logger.warning("PAZy reported count (%r) does not match number of entries read (%d)", rawD["count"], len(entryL))
        else:
            entryL = rawD
        #
        oD = {}
        for entry in entryL:
            try:
                enzymeId = str(entry["id"]).strip()
                sequence = self.__cleanSequence(entry.get("amino_acid_sequence"))
                if not sequence:
                    logger.debug("Skipping PAZy enzyme %r without sequence", enzymeId)
                    continue
                substrateL = []
                for sD in entry.get("substrates") or []:
                    sName = (sD.get("abbreviation") or sD.get("name") or "").strip()
                    if sName and sName not in substrateL:
                        substrateL.append(sName)
                eD = {
                    "source_db": "PAZy",
                    "name": (entry.get("name") or "").strip(),
                    "substrates": sorted(substrateL),
                    "sequence": sequence,
                }
                orgD = entry.get("organism") or {}
                orgName = (orgD.get("scientific_name") or "").strip()
                if orgName:
                    eD["organism_scientific_name"] = orgName
                taxId = orgD.get("ncbi_taxonomy_id")
                if taxId is not None and str(taxId).strip():
                    eD["organism_ncbi_tax_id"] = str(taxId).strip()
                uniprotId = (entry.get("uniprot_accession") or "").strip()
                if uniprotId:
                    eD["uniprot_accession"] = uniprotId
                oD[enzymeId] = eD
            except Exception as e:
                logger.exception("Failing to process PAZy entry %r with %s", entry.get("id") if isinstance(entry, dict) else entry, str(e))
        #
        if len(oD) < minCount:
            self.__fU.remove(downloadPath)
            logger.error("Too few PAZy records (%d < %d) in the fetched data", len(oD), minCount)
            return None, None
        return oD, downloadPath

    def __fetchPlasticDBData(self, minCount=100):
        """Fetch PlasticDB FASTA data and restructure it (written to plasticdb-select-data.json by fetchPlasticEnzymeData())

        PlasticDB fasta headers have the form: >{enzymeId}|{name}|{organism}|{plastics}|{accession}
        (e.g., >00001|PHB_depolymerase|Ralstonia_pickettii|PHB_PHA|BAA04986.1), where spaces in the name and organism
        are written as "_", plastics are "_" separated (spaces within a plastic name are written as "-", e.g. "PBS-Blend"),
        and a missing accession is "NA".

        Returns:
            (dict, str): standardized records keyed by enzyme ID and the path to the downloaded raw file, or (None, None) if
                         the data couldn't be fetched or has too few records (the download is then removed)
        """
        downloadPath = self.__getSourceFile("PlasticDB", self.__plasticDBUrl)
        if not downloadPath:
            return None, None
        # The "default" comment style returns {header (without ">"): {"sequence": ...}}
        sD = self.__mU.doImport(downloadPath, fmt="fasta", commentStyle="default")
        if not sD:
            self.__fU.remove(downloadPath)
            logger.error("Unable to read PlasticDB data from %r", self.__plasticDBUrl)
            return None, None
        #
        oD = {}
        for header, seqD in sD.items():
            try:
                fieldL = [f.strip() for f in header.split("|")]
                if len(fieldL) != 5:
                    logger.warning("Skipping PlasticDB record with unexpected header format %r", header)
                    continue
                enzymeId, name, organism, plastics, accession = fieldL
                sequence = self.__cleanSequence(seqD.get("sequence"))
                if not enzymeId or not sequence:
                    logger.debug("Skipping PlasticDB record %r without enzyme ID or sequence", header)
                    continue
                if enzymeId in oD:
                    logger.warning("Duplicate PlasticDB enzyme ID %r - keeping the first", enzymeId)
                    continue
                substrateL = []
                for substrate in plastics.split("_"):
                    if substrate and substrate not in substrateL:
                        substrateL.append(substrate)
                eD = {
                    "source_db": "PlasticDB",
                    "name": name.replace("_", " ").strip(),
                    "substrates": sorted(substrateL),
                    "sequence": sequence,
                }
                orgName = organism.replace("_", " ").strip()
                if orgName:
                    eD["organism_scientific_name"] = orgName
                # PlasticDB doesn't provide taxonomy IDs in the FASTA file. Accessions are mostly GenBank protein IDs,
                # with some UniProt accessions (e.g., Q96UT0.1)
                if accession and accession.upper() != "NA":
                    if self.__isUniProtAccession(accession.split(".")[0]):
                        eD["uniprot_accession"] = accession.split(".")[0]
                    else:
                        eD["genbank_accession"] = accession
                oD[enzymeId] = eD
            except Exception as e:
                logger.exception("Failing to process PlasticDB record %r with %s", header, str(e))
        #
        if len(oD) < minCount:
            self.__fU.remove(downloadPath)
            logger.error("Too few PlasticDB records (%d < %d) in the fetched data", len(oD), minCount)
            return None, None
        return oD, downloadPath

    def __isUniProtAccession(self, accession):
        return re.match(r"^([OPQ][0-9][A-Z0-9]{3}[0-9]|[A-NR-Z][0-9]([A-Z][A-Z0-9]{2}[0-9]){1,2})$", accession) is not None

    def __cleanSequence(self, sequence):
        """Remove all whitespace from a raw sequence string (returns an empty string if no sequence)."""
        if not sequence:
            return ""
        return re.sub(r"\s+", "", str(sequence)).upper()

    def __exportSelectData(self, filePath, dataD, sourceDb):
        tS = datetime.datetime.now().isoformat()
        # Version is hard-coded since it is not available via the source databases
        sD = {"version": "1.0", "created": tS, "data": dataD}
        ok = self.__mU.doExport(filePath, sD, fmt="json", indent=3)
        logger.info("Export %s data (%d) status %r", sourceDb, len(dataD), ok)
        return sD if ok else {}
