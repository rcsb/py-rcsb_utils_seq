##
#  File:           PlasticEnzymeAnnotationProvider.py
#  Date:           7-Oct-2026
#
#  Updates:
#
##
"""
Accessors and generators for plastic-degrading enzyme annotation data (sourced from PAZy and PlasticDB).
"""

import datetime
import logging
import os.path
import re
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
        #
        # Source locators (remote URLs by default; local file paths are also accepted by FileUtil.get())
        self.__pazyUrl = kwargs.get("pazyUrl", "https://api.pazy.eu/api/proteins/json/?format=json")
        self.__plasticDBUrl = kwargs.get("plasticDBUrl", "https://plasticdb.org/static/degraders_list.tsv")
        #
        self.__mU = MarshalUtil(workPath=self.__dirPath)
        self.__fU = FileUtil()
        self.__fU.mkdir(self.__dirPath)

    def __getPazyDataPath(self):
        return os.path.join(self.__dirPath, "pazy-select-data.json")

    def __getPlasticDBDataPath(self):
        return os.path.join(self.__dirPath, "plasticdb-select-data.json")

    def fetchPlasticEnzymeData(self):
        """Fetch the PAZy and PlasticDB source data and restructure each into a standardized JSON file
        (pazy-select-data.json and plasticdb-select-data.json) in the cache path.

        Returns:
            bool: True for success (both sources) or False otherwise
        """
        startTime = time.time()
        pazyOk = self.__fetchPazyData()
        plasticDBOk = self.__fetchPlasticDBData()
        ok = pazyOk and plasticDBOk
        logger.info(
            "Completed fetch of plastic enzyme data (%r) PAZy (%r) PlasticDB (%r) at %s (%.4f seconds)",
            ok,
            pazyOk,
            plasticDBOk,
            time.strftime("%Y %m %d %H:%M:%S", time.localtime()),
            time.time() - startTime,
        )
        return ok

    def exportFasta(self, fastaPath=None):
        """Export a combined (PAZy and PlasticDB) plastic enzyme sequence fasta file (plastic-enzyme-sequences.fa)
        using the standardized select-data JSON files in the cache path.

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
            sD = {}
            for dataPath in [self.__getPazyDataPath(), self.__getPlasticDBDataPath()]:
                if not self.__mU.exists(dataPath):
                    logger.error("Missing plastic enzyme data file %r (run fetchPlasticEnzymeData first)", dataPath)
                    return False
                qD = self.__mU.doImport(dataPath, fmt="json")
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

    def __getFastaPath(self):
        return os.path.join(self.__dirPath, "plastic-enzyme-sequences.fa")

    def exportFastaFromFastaFiles(self, pazyFastaPath="pazy_proteins.fasta", plasticDBFastaPath="PlasticDB.fasta", fastaPath=None):
        """Export a combined plastic enzyme sequence fasta file (fa-plastic-enzymes-sequences.fa) using only the
        fasta files downloaded from PAZy and PlasticDB (i.e., not the JSON/TSV data files).

        The output has the same format as exportFasta(): >{enzymeId}|enzymeId|{sourceDb}|sourceDb

        PAZy fasta headers have the form ">{enzymeId} {name}" (e.g., ">1 IsPETase"), and
        PlasticDB fasta headers have the form ">{enzymeId}|{name}|{organism}|{plastics}|{accession}" (e.g., ">00001|PHB_depolymerase|...").

        Args:
            pazyFastaPath (str): PAZy fasta input file path
            plasticDBFastaPath (str): PlasticDB fasta input file path
            fastaPath (str, optional): fasta output file path. Defaults to fa-plastic-enzymes-sequences.fa in the cache path

        Returns:
            bool: True for success or False otherwise
        """
        fastaPath = fastaPath if fastaPath else os.path.join(self.__dirPath, "fa-plastic-enzymes-sequences.fa")
        pazyDataPath = os.path.join(self.__dirPath, "pazy-fasta-select-data.json")
        plasticDBDataPath = os.path.join(self.__dirPath, "plasticdb-fasta-select-data.json")
        ok = False
        try:
            # Step 1: restructure each source fasta file into a standardized select-data JSON file
            pazyOk = self.__buildSelectDataFromFasta(pazyFastaPath, "PAZy", pazyDataPath)
            plasticDBOk = self.__buildSelectDataFromFasta(plasticDBFastaPath, "PlasticDB", plasticDBDataPath)
            if not (pazyOk and plasticDBOk):
                return False
            # Step 2: create the combined fasta file from the select-data JSON files
            sD = {}
            for dataPath in [pazyDataPath, plasticDBDataPath]:
                qD = self.__mU.doImport(dataPath, fmt="json")
                for enzymeId, tD in qD["data"].items():
                    if not tD.get("sequence"):
                        continue
                    # Dictionary order determines the fasta header: {value}|{key}|{value}|{key}
                    cId = "|".join([str(enzymeId), "enzymeId", str(tD["source_db"]), "sourceDb"])
                    sD[cId] = {"sequence": tD["sequence"], "enzymeId": enzymeId, "sourceDb": tD["source_db"]}
            #
            ok = self.__mU.doExport(fastaPath, sD, fmt="fasta", makeComment=True)
            logger.info("Export plastic enzyme fasta from fasta files (%d) status %r", len(sD), ok)
        except Exception as e:
            logger.exception("Failing to export plastic enzyme fasta from fasta files with %s", str(e))
        return ok

    def __buildSelectDataFromFasta(self, inPath, sourceDb, dataPath):
        """Restructure a PAZy or PlasticDB fasta file into a standardized select-data JSON file. Only the enzyme ID,
        source database, name and sequence are available from the fasta files.

        Args:
            inPath (str): source fasta file path
            sourceDb (str): "PAZy" or "PlasticDB"
            dataPath (str): select-data JSON output file path

        Returns:
            bool: True for success or False otherwise
        """
        if not self.__mU.exists(inPath):
            logger.error("Missing %s fasta file %r", sourceDb, inPath)
            return False
        oD = {}
        for header, sequence in self.__readFastaRecords(inPath):
            if sourceDb == "PlasticDB":
                # PlasticDB headers are "|" separated: {enzymeId}|{name}|{organism}|{plastics}|{accession}
                fieldL = [f.strip() for f in header.split("|")]
                enzymeId = fieldL[0]
                name = fieldL[1] if len(fieldL) > 1 else ""
            else:
                # PAZy headers are: {enzymeId} {name}
                fieldL = header.split(None, 1)
                enzymeId = fieldL[0] if fieldL else ""
                name = fieldL[1].strip() if len(fieldL) > 1 else ""
            if not enzymeId or not sequence:
                logger.debug("Skipping %s fasta record %r without enzyme ID or sequence", sourceDb, header)
                continue
            if enzymeId in oD:
                logger.warning("Duplicate %s enzyme ID %r in fasta file - keeping the first", sourceDb, enzymeId)
                continue
            oD[enzymeId] = {"source_db": sourceDb, "name": name, "sequence": sequence}
        return self.__exportSelectData(dataPath, oD, sourceDb)

    def __readFastaRecords(self, filePath):
        """Read a fasta file and return a list of (header, sequence) tuples, where the header excludes the leading ">"
        and the sequence has all whitespace removed and is upper case.
        """
        recordL = []
        header = None
        seqL = []
        with open(filePath, "r", encoding="utf-8") as ifh:
            for line in ifh:
                line = line.rstrip("\r\n")
                if line.startswith(">"):
                    if header is not None:
                        recordL.append((header, self.__cleanSequence("".join(seqL))))
                    header, seqL = line[1:].strip(), []
                elif header is not None:
                    seqL.append(line)
        if header is not None:
            recordL.append((header, self.__cleanSequence("".join(seqL))))
        return recordL

    def __fetchPazyData(self):
        """Fetch PAZy JSON data and restructure it -> pazy-select-data.json

        Returns:
            bool: True for success or False otherwise
        """
        rawPath = os.path.join(self.__dirPath, "pazy_proteins.json")
        logger.info("Fetching PAZy data from %r to %r", self.__pazyUrl, rawPath)
        if not self.__fU.get(self.__pazyUrl, rawPath):
            logger.error("Fetching failed for PAZy data")
            return False
        rawD = self.__mU.doImport(rawPath, fmt="json")
        if not rawD:
            logger.error("Unable to read PAZy data from %r", rawPath)
            return False
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
                    "substrates": substrateL,
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
        return self.__exportSelectData(self.__getPazyDataPath(), oD, "PAZy")

    def __fetchPlasticDBData(self):
        """Fetch PlasticDB TSV data and restructure it -> plasticdb-select-data.json

        Multiple rows may describe the same enzyme (one row per plastic), so rows are combined
        into a single record per enzyme ID with all plastics collected in the substrates list.

        Returns:
            bool: True for success or False otherwise
        """
        rawPath = os.path.join(self.__dirPath, "degraders_list.tsv")
        logger.info("Fetching PlasticDB data from %r to %r", self.__plasticDBUrl, rawPath)
        if not self.__fU.get(self.__plasticDBUrl, rawPath):
            logger.error("Fetching failed for PlasticDB data")
            return False
        # uncomment=False so that "#" characters and empty trailing columns in the data are preserved
        rowL = self.__mU.doImport(rawPath, fmt="tdd", rowFormat="dict", uncomment=False)
        if not rowL:
            logger.error("Unable to read PlasticDB data from %r", rawPath)
            return False
        #
        oD = {}
        for rowD in rowL:
            try:
                # Note that the source header contains a typo ("Enzyne")
                dbName = rowD.get("Database Enzyne Name", rowD.get("Database Enzyme Name", "")) or ""
                enzymeId = dbName.split("|")[0].strip()
                sequence = self.__cleanSequence(rowD.get("Sequence"))
                if not enzymeId or not sequence:
                    continue
                substrate = (rowD.get("Plastic") or "").strip()
                if enzymeId in oD:
                    eD = oD[enzymeId]
                    if eD["sequence"] != sequence:
                        logger.warning("PlasticDB enzyme %r has multiple different sequences - keeping the first", enzymeId)
                        continue
                    if substrate and substrate not in eD["substrates"]:
                        eD["substrates"].append(substrate)
                    continue
                eD = {
                    "source_db": "PlasticDB",
                    "name": (rowD.get("Enzyme") or "").strip(),
                    "substrates": [substrate] if substrate else [],
                    "sequence": sequence,
                }
                orgName = (rowD.get("Microorganism") or "").strip()
                if orgName:
                    eD["organism_scientific_name"] = orgName
                taxId = (rowD.get("Tax ID") or "").strip()
                if taxId:
                    eD["organism_ncbi_tax_id"] = taxId
                # PlasticDB provides GenBank (not UniProt) accessions, so "uniprot_accession" is not available
                oD[enzymeId] = eD
            except Exception as e:
                logger.exception("Failing to process PlasticDB row %r with %s", rowD, str(e))
        #
        return self.__exportSelectData(self.__getPlasticDBDataPath(), oD, "PlasticDB")

    def __cleanSequence(self, sequence):
        """Remove all whitespace from a raw sequence string (returns an empty string if no sequence)."""
        if not sequence:
            return ""
        return re.sub(r"\s+", "", str(sequence)).upper()

    def __exportSelectData(self, filePath, dataD, sourceDb):
        tS = datetime.datetime.now().isoformat()
        # Version is hard-coded since it is not available via the source databases
        ok = self.__mU.doExport(filePath, {"version": "1.0", "created": tS, "data": dataD}, fmt="json", indent=3)
        logger.info("Export %s data (%d) status %r", sourceDb, len(dataD), ok)
        return ok
