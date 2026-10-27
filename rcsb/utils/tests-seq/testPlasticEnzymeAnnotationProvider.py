##
# File:    testPlasticEnzymeAnnotationProvider.py
# Date:    20-Oct-2026
#
# Updates:
#  26-Oct-2026 dwp Add testReloadMalformedCache
#  27-Oct-2026 dwp Replace the opt-in live fetch test with tests that always check that the PAZy and PlasticDB
#                  source data are accessible and usable (each fails if its source can't be fetched); start each test
#                  from an empty cache; add tests for the cached-data and stash fallbacks, fetching with useCache,
#                  and the TEMP local raw source files
##
"""
Tests for the plastic-degrading enzyme (PAZy and PlasticDB) annotation provider.
"""

import json
import logging
import os
import shutil
import unittest
from unittest import mock

from rcsb.utils.io.MarshalUtil import MarshalUtil
from rcsb.utils.seq.PlasticEnzymeAnnotationProvider import PlasticEnzymeAnnotationProvider

HERE = os.path.abspath(os.path.dirname(__file__))
TOPDIR = os.path.dirname(os.path.dirname(os.path.dirname(HERE)))

logging.basicConfig(level=logging.INFO, format="%(asctime)s [%(levelname)s]-%(module)s.%(funcName)s: %(message)s")
logger = logging.getLogger()

SEQ_A = "MRNTLKAAFKLGVISAALLAPFATQAATAGPGAWSSQQTWAADSVNGG"
SEQ_B = "MKHPYGYRWHWLYALVVTLMTALATFSAHAAVTAGPGAWSSQQTWAADTV"

PAZY_JSON = {
    "count": 3,
    "results": [
        {
            "id": 15,
            "name": "RgPETase",
            "amino_acid_sequence": SEQ_A[:24] + "\n" + SEQ_A[24:],
            "organism": {"scientific_name": "Rhizobacter gummiphilus", "ncbi_taxonomy_id": 946333},
            "uniprot_accession": "A0A1W6L588",
            "substrates": [{"abbreviation": "PET", "name": "Polyethylene terephthalate"}, {"abbreviation": "PET", "name": "Polyethylene terephthalate"}],
        },
        {"id": 16, "name": "NoSequence", "amino_acid_sequence": None, "organism": None, "uniprot_accession": None, "substrates": []},
        {"id": 17, "name": "PHB depolymerase", "amino_acid_sequence": SEQ_B, "organism": {}, "uniprot_accession": None, "substrates": [{"abbreviation": "PHA"}]},
    ],
}

PLASTICDB_FASTA = (
    ">00001|PHB_depolymerase|Ralstonia_pickettii|PHB_PHA|BAA04986.1\n"
    + SEQ_B
    + "\n>00008|Cutinase|Fusarium_verticillioides|PCL|Q96UT0.1\n"
    + SEQ_A
    + "\n>00068|PU_esterase_|Pseudozyma_antarctica|PBS-Blend_O-PVA|NA\n"
    + SEQ_A[::-1]
    + "\n>00069|BadHeader|NoPlastics\n"
    + SEQ_B
    + "\n"
)


class PlasticEnzymeAnnotationProviderTests(unittest.TestCase):
    def setUp(self):
        self.__workPath = os.path.join(HERE, "test-output", "plastic-enzyme")
        self.__cachePath = os.path.join(self.__workPath, "CACHE")
        shutil.rmtree(self.__cachePath, ignore_errors=True)
        self.__mU = MarshalUtil(workPath=self.__workPath)
        self.__mU.mkdir(self.__workPath)
        self.__pazyPath = os.path.join(self.__workPath, "pazy-test.json")
        self.__plasticDBPath = os.path.join(self.__workPath, "plasticdb-test.fasta")
        self.__mU.doExport(self.__pazyPath, PAZY_JSON, fmt="json")
        with open(self.__plasticDBPath, "w", encoding="utf-8") as ofh:
            ofh.write(PLASTICDB_FASTA)

    def tearDown(self):
        pass

    def __getTestProvider(self, useCache=False):
        return PlasticEnzymeAnnotationProvider(cachePath=self.__cachePath, useCache=useCache, pazyUrl=self.__pazyPath, plasticDBUrl=self.__plasticDBPath)

    def testFetchAndStandardize(self):
        """PAZy and PlasticDB retrieval and normalization (local test files)."""
        peP = self.__getTestProvider()
        self.assertFalse(peP.testSourceCache(minCount=0))
        ok = peP.fetchPlasticEnzymeData(minCount=1)
        self.assertTrue(ok)
        self.assertTrue(peP.testSourceCache(minCount=1))
        #
        # PAZy: entries without a sequence are skipped; sequences are cleaned; substrates are unique abbreviations
        self.assertIsNone(peP.getEnzyme("PAZy", "16"))
        pD = peP.getEnzyme("PAZy", "15")
        self.assertEqual(pD["sequence"], SEQ_A)
        self.assertEqual(pD["substrates"], ["PET"])
        self.assertEqual(pD["organism_ncbi_tax_id"], "946333")
        self.assertEqual(pD["uniprot_accession"], "A0A1W6L588")
        self.assertNotIn("organism_scientific_name", peP.getEnzyme("PAZy", "17"))
        #
        # PlasticDB: all header fields parsed; records with unexpected headers are skipped
        self.assertIsNone(peP.getEnzyme("PlasticDB", "00069"))
        # Downloaded source files are removed after parsing (and never left as a PlasticDB.fasta for the TEMP local-file check)
        self.assertEqual(sorted(os.listdir(os.path.join(self.__cachePath, "plastic-enzyme-annotations"))), ["pazy-select-data.json", "plasticdb-select-data.json"])
        pD = peP.getEnzyme("PlasticDB", "00001")
        self.assertEqual(pD["name"], "PHB depolymerase")
        self.assertEqual(pD["organism_scientific_name"], "Ralstonia pickettii")
        self.assertEqual(pD["substrates"], ["PHA", "PHB"])
        self.assertEqual(pD["genbank_accession"], "BAA04986.1")
        self.assertEqual(peP.getEnzyme("PlasticDB", "00008")["uniprot_accession"], "Q96UT0")
        pD = peP.getEnzyme("PlasticDB", "00068")
        self.assertEqual(pD["name"], "PU esterase")
        self.assertEqual(pD["substrates"], ["O-PVA", "PBS-Blend"])
        self.assertNotIn("genbank_accession", pD)
        self.assertNotIn("uniprot_accession", pD)

    def testFetchMinCount(self):
        """Source data with too few records doesn't replace the cache."""
        peP = self.__getTestProvider()
        ok = peP.fetchPlasticEnzymeData(minCount=100)
        self.assertFalse(ok)
        self.assertFalse(peP.testSourceCache(minCount=0))

    def testFetchFailure(self):
        """A failed source fetch is logged as an error and doesn't write the source data."""
        peP = PlasticEnzymeAnnotationProvider(cachePath=self.__cachePath, useCache=False, pazyUrl=os.path.join(self.__workPath, "missing.json"), plasticDBUrl=self.__plasticDBPath)
        with self.assertLogs("rcsb.utils.seq.PlasticEnzymeAnnotationProvider", level="ERROR") as cm:
            ok = peP.fetchPlasticEnzymeData(minCount=1)
        self.assertFalse(ok)
        self.assertTrue(any("Fetching PAZy data from" in msg for msg in cm.output))
        self.assertIsNone(peP.getEnzyme("PAZy", "15"))
        self.assertIsNotNone(peP.getEnzyme("PlasticDB", "00001"))
        self.assertTrue(any("no stash configuration given" in msg for msg in cm.output))
        self.assertFalse([fn for fn in os.listdir(os.path.join(self.__cachePath, "plastic-enzyme-annotations")) if fn.endswith(".download")])

    def testFetchFallbackToCache(self):
        """A source that can't be fetched keeps its existing cached data (without a stash restore)."""
        self.assertTrue(self.__getTestProvider().fetchPlasticEnzymeData(minCount=1))
        peP = PlasticEnzymeAnnotationProvider(cachePath=self.__cachePath, useCache=False, pazyUrl=os.path.join(self.__workPath, "missing.json"), plasticDBUrl=self.__plasticDBPath)
        with mock.patch.object(PlasticEnzymeAnnotationProvider, "restore") as mRestore:
            ok = peP.fetchPlasticEnzymeData(minCount=1, cfgOb="cfgOb", configName="configName")
        self.assertFalse(ok)
        mRestore.assert_not_called()
        self.assertEqual(peP.getEnzyme("PAZy", "15")["sequence"], SEQ_A)
        self.assertIsNotNone(peP.getEnzyme("PlasticDB", "00001"))
        self.assertTrue(peP.testSourceCache(minCount=1))

    def testFetchFallbackToStash(self):
        """A source that can't be fetched or found in the cache is restored from the stash, without replacing the freshly fetched source."""
        dirPath = os.path.join(self.__cachePath, "plastic-enzyme-annotations")

        def fakeRestore(cfgOb, configName, remotePrefix=None, useStash=True, useGit=False):
            # The stashed backup holds older data for both sources
            self.assertEqual((cfgOb, configName, remotePrefix, useStash, useGit), ("cfgOb", "configName", "A", True, False))
            for fileName, sourceDb, enzymeId in [("pazy-select-data.json", "PAZy", "99"), ("plasticdb-select-data.json", "PlasticDB", "99999")]:
                eD = {"source_db": sourceDb, "name": "stashed", "substrates": ["PET"], "sequence": SEQ_A}
                self.__mU.doExport(os.path.join(dirPath, fileName), {"version": "1.0", "created": "stash", "data": {enzymeId: eD, enzymeId + "0": eD}}, fmt="json")
            return True

        peP = PlasticEnzymeAnnotationProvider(cachePath=self.__cachePath, useCache=False, pazyUrl=os.path.join(self.__workPath, "missing.json"), plasticDBUrl=self.__plasticDBPath)
        with mock.patch.object(PlasticEnzymeAnnotationProvider, "restore", side_effect=fakeRestore) as mRestore:
            ok = peP.fetchPlasticEnzymeData(minCount=1, cfgOb="cfgOb", configName="configName", remotePrefix="A")
        self.assertFalse(ok)
        mRestore.assert_called_once()
        # PAZy comes from the stash; PlasticDB is the freshly fetched data (not the stashed copy)
        self.assertEqual(peP.getEnzyme("PAZy", "99")["name"], "stashed")
        self.assertIsNotNone(peP.getEnzyme("PlasticDB", "00001"))
        self.assertIsNone(peP.getEnzyme("PlasticDB", "99999"))
        self.assertTrue(peP.testSourceCache(minCount=1))

    def testFetchUseCache(self):
        """With useCache, valid cached source data is used without fetching; without it, the sources are fetched."""
        self.assertTrue(self.__getTestProvider().fetchPlasticEnzymeData(minCount=1))
        missingKw = {"pazyUrl": os.path.join(self.__workPath, "missing.json"), "plasticDBUrl": os.path.join(self.__workPath, "missing.fasta")}
        peP = PlasticEnzymeAnnotationProvider(cachePath=self.__cachePath, useCache=False, **missingKw)
        with self.assertNoLogs("rcsb.utils.seq.PlasticEnzymeAnnotationProvider", level="ERROR"):
            self.assertTrue(peP.fetchPlasticEnzymeData(minCount=1, useCache=True))
        self.assertIsNotNone(peP.getEnzyme("PAZy", "15"))
        self.assertIsNotNone(peP.getEnzyme("PlasticDB", "00001"))
        # Cached data with too few records isn't used (and here can't be replaced, since the sources are missing)
        with self.assertLogs("rcsb.utils.seq.PlasticEnzymeAnnotationProvider", level="ERROR"):
            self.assertFalse(peP.fetchPlasticEnzymeData(minCount=100, useCache=True))
        # Without useCache the sources are always fetched (here failing, so the cached data is kept)
        with self.assertLogs("rcsb.utils.seq.PlasticEnzymeAnnotationProvider", level="ERROR"):
            self.assertFalse(peP.fetchPlasticEnzymeData(minCount=1))
        self.assertTrue(peP.testSourceCache(minCount=1))

    def testTempLocalSourceFiles(self):
        """TEMP: with useCache, raw source files already in the cache directory are used instead of fetching, and kept."""
        dirPath = os.path.join(self.__cachePath, "plastic-enzyme-annotations")
        self.__mU.mkdir(dirPath)
        shutil.copyfile(self.__pazyPath, os.path.join(dirPath, "pazy_proteins.json"))
        shutil.copyfile(self.__plasticDBPath, os.path.join(dirPath, "PlasticDB.fasta"))
        missingKw = {"pazyUrl": os.path.join(self.__workPath, "missing.json"), "plasticDBUrl": os.path.join(self.__workPath, "missing.fasta")}
        # Without useCache, the raw files are ignored and the sources are fetched (here failing)
        peP = PlasticEnzymeAnnotationProvider(cachePath=self.__cachePath, useCache=False, **missingKw)
        self.assertFalse(peP.fetchPlasticEnzymeData(minCount=1))
        self.assertFalse(peP.testSourceCache(minCount=0))
        # With useCache, both raw files are used and left in place
        self.assertTrue(peP.fetchPlasticEnzymeData(minCount=1, useCache=True))
        self.assertEqual(peP.getEnzyme("PAZy", "15")["sequence"], SEQ_A)
        self.assertEqual(peP.getEnzyme("PlasticDB", "00001")["name"], "PHB depolymerase")
        self.assertTrue(os.path.exists(os.path.join(dirPath, "pazy_proteins.json")))
        self.assertTrue(os.path.exists(os.path.join(dirPath, "PlasticDB.fasta")))

    def testFetchStashKeepsCachedData(self):
        """A stash restore for one source keeps the valid cached data in use for the other source."""
        dirPath = os.path.join(self.__cachePath, "plastic-enzyme-annotations")
        self.assertTrue(self.__getTestProvider().fetchPlasticEnzymeData(minCount=1))
        os.remove(os.path.join(dirPath, "plasticdb-select-data.json"))

        def fakeRestore(cfgOb, configName, remotePrefix=None, useStash=True, useGit=False):
            _ = (cfgOb, configName, remotePrefix, useStash, useGit)
            for fileName, sourceDb, enzymeId in [("pazy-select-data.json", "PAZy", "99"), ("plasticdb-select-data.json", "PlasticDB", "99999")]:
                eD = {"source_db": sourceDb, "name": "stashed", "substrates": ["PET"], "sequence": SEQ_A}
                self.__mU.doExport(os.path.join(dirPath, fileName), {"version": "1.0", "created": "stash", "data": {enzymeId: eD, enzymeId + "0": eD}}, fmt="json")
            return True

        peP = PlasticEnzymeAnnotationProvider(cachePath=self.__cachePath, useCache=False, pazyUrl=self.__pazyPath, plasticDBUrl=os.path.join(self.__workPath, "missing.fasta"))
        with mock.patch.object(PlasticEnzymeAnnotationProvider, "restore", side_effect=fakeRestore) as mRestore:
            self.assertFalse(peP.fetchPlasticEnzymeData(minCount=1, useCache=True, cfgOb="cfgOb", configName="configName"))
        mRestore.assert_called_once()
        # PAZy keeps its cached data (not the stashed copy); PlasticDB comes from the stash
        self.assertIsNotNone(peP.getEnzyme("PAZy", "15"))
        self.assertIsNone(peP.getEnzyme("PAZy", "99"))
        self.assertEqual(peP.getEnzyme("PlasticDB", "99999")["name"], "stashed")

    def testExportFasta(self):
        """FASTA generation and header format."""
        peP = self.__getTestProvider()
        self.assertFalse(peP.exportFasta())
        self.assertTrue(peP.fetchPlasticEnzymeData(minCount=1))
        fastaPath = os.path.join(self.__workPath, "plastic-enzyme-sequences-test.fa")
        self.assertTrue(peP.exportFasta(fastaPath))
        sD = self.__mU.doImport(fastaPath, fmt="fasta", commentStyle="default")
        expectedL = ["15|PAZy", "17|PAZy", "00001|PlasticDB", "00008|PlasticDB", "00068|PlasticDB"]
        self.assertEqual(sorted(sD.keys()), sorted(["%s|enzymeId|%s|sourceDb" % tuple(ky.split("|")) for ky in expectedL]))
        self.assertEqual(sD["15|enzymeId|PAZy|sourceDb"]["sequence"], SEQ_A)

    def testBuildAnnotationList(self):
        """Parsing of sequence search results, mapping back to source annotations, cache reload and entity retrieval."""
        peP = self.__getTestProvider()
        self.assertFalse(peP.buildAnnotationList(os.path.join(self.__workPath, "missing.json"), minCount=1))
        self.assertTrue(peP.fetchPlasticEnzymeData(minCount=1))
        #
        tgt = "{}|entityId|1|srcId|1|seqBeg|292|seqEnd|292|seqLen|9606|taxId"
        mD = {
            "15|enzymeId|PAZy|sourceDb": [
                {"query": "15|enzymeId|PAZy|sourceDb", "target": tgt.format("7DZT_1"), "sequenceIdentity": 100.0, "bitScore": 752.0},
                {"query": "15|enzymeId|PAZy|sourceDb", "target": tgt.format("pdb_00001abc_2"), "sequenceIdentity": 97.5, "bitScore": 700.0},
            ],
            "00008|enzymeId|PlasticDB|sourceDb": [
                {"query": "00008|enzymeId|PlasticDB|sourceDb", "target": tgt.format("7DZT_1"), "sequenceIdentity": 100.0, "bitScore": 752.0},
                # Second source part of the same entity: only the best match per (entity, enzyme) is kept
                {"query": "00008|enzymeId|PlasticDB|sourceDb", "target": tgt.format("7DZT_1").replace("1|srcId", "2|srcId"), "sequenceIdentity": 96.0, "bitScore": 300.0},
            ],
            "00001|enzymeId|PlasticDB|sourceDb": [
                {"query": "00001|enzymeId|PlasticDB|sourceDb", "target": tgt.format("7DZT_1"), "sequenceIdentity": 99.0, "bitScore": 750.0},
            ],
            "99999|enzymeId|PlasticDB|sourceDb": [
                {"query": "99999|enzymeId|PlasticDB|sourceDb", "target": tgt.format("9XYZ_1"), "sequenceIdentity": 100.0, "bitScore": 10.0},
            ],
        }
        matchPath = os.path.join(self.__workPath, "plastic-enzyme-vs-pdbprent-filtered-results-test.json")
        with open(matchPath, "w", encoding="utf-8") as ofh:
            json.dump(mD, ofh)
        self.assertTrue(peP.buildAnnotationList(matchPath, minCount=1))
        #
        aDL = peP.getAnnotation("7DZT_1")
        self.assertEqual([(aD["source_db"], aD["annotation_id"]) for aD in aDL], [("PAZy", "15"), ("PlasticDB", "00008"), ("PlasticDB", "00001")])
        self.assertEqual([aD["perfect_match"] for aD in aDL], ["Y", "Y", "N"])
        aD = aDL[0]
        self.assertEqual((aD["entry_id"], aD["entity_id"], aD["name"], aD["substrates"]), ("7DZT", "1", "RgPETase", ["PET"]))
        self.assertEqual(aD["description"], "Plastic degrading enzyme (sequence matched against PAZy).")
        self.assertEqual(aD["reference_scheme"], "PDB entity")
        self.assertEqual(aDL[1]["seq_id_pct"], 100.0)
        #
        aDL = peP.getAnnotation("pdb_00001abc_2")
        self.assertEqual((aDL[0]["entry_id"], aDL[0]["entity_id"]), ("pdb_00001abc", "2"))
        self.assertIsNone(peP.getAnnotation("9XYZ_1"))
        self.assertFalse(peP.hasAnnotation("7dzt_1"))
        self.assertTrue(peP.testCache(minCount=1))
        self.assertFalse(peP.testCache(minCount=100))
        #
        # Empty or missing search results don't overwrite the existing annotation data
        emptyPath = os.path.join(self.__workPath, "empty-results.json")
        with open(emptyPath, "w", encoding="utf-8") as ofh:
            json.dump({}, ofh)
        self.assertFalse(peP.buildAnnotationList(emptyPath, minCount=1))
        self.assertFalse(peP.buildAnnotationList(os.path.join(self.__workPath, "missing.json"), minCount=1))
        self.assertFalse(peP.buildAnnotationList(matchPath, minCount=100))
        #
        # Cache reload
        peP = self.__getTestProvider(useCache=True)
        self.assertTrue(peP.testCache(minCount=1))
        self.assertTrue(peP.hasAnnotation("7DZT_1"))
        self.assertEqual(peP.getVersion(), "1.0")
        # useCache=False loads nothing
        peP = self.__getTestProvider(useCache=False)
        self.assertFalse(peP.testCache(minCount=0))
        self.assertIsNone(peP.getAnnotation("7DZT_1"))

    def testReloadMalformedCache(self):
        """Cache files that don't load as a dictionary are ignored."""
        cachePath = os.path.join(self.__workPath, "CACHE-malformed")
        dirPath = os.path.join(cachePath, "plastic-enzyme-annotations")
        self.__mU.mkdir(dirPath)
        for fileName in ["pazy-select-data.json", "plasticdb-select-data.json", "plastic-enzyme-annotation-data.json"]:
            self.__mU.doExport(os.path.join(dirPath, fileName), ["not", "a", "dictionary"], fmt="json")
        with self.assertLogs("rcsb.utils.seq.PlasticEnzymeAnnotationProvider", level="WARNING") as cm:
            peP = PlasticEnzymeAnnotationProvider(cachePath=cachePath, useCache=True)
        self.assertEqual(sum("Unable to load plastic enzyme cache file" in msg for msg in cm.output), 3)
        self.assertFalse(peP.testCache(minCount=0))
        self.assertIsNone(peP.getVersion())
        self.assertIsNone(peP.getAnnotation("7DZT_1"))
        self.assertIsNone(peP.getEnzyme("PAZy", "15"))

    def __checkLiveSource(self, sourceDb, urlKey, dataFileName, minCount=100):
        """Fetch one source from its default (live) URL, and the other from its local test file, and check that the live
        source returns enough usable records."""
        cachePath = os.path.join(self.__workPath, "CACHE-live-%s" % sourceDb)
        shutil.rmtree(cachePath, ignore_errors=True)
        kwargs = {"pazyUrl": self.__pazyPath, "plasticDBUrl": self.__plasticDBPath}
        kwargs.pop(urlKey)
        peP = PlasticEnzymeAnnotationProvider(cachePath=cachePath, useCache=False, **kwargs)
        ok = peP.fetchPlasticEnzymeData(minCount=1)
        self.assertTrue(ok, "%s source data could not be fetched (see the error logged above)" % sourceDb)
        dataD = self.__mU.doImport(os.path.join(cachePath, "plastic-enzyme-annotations", dataFileName), fmt="json")["data"]
        self.assertGreater(len(dataD), minCount, "too few %s records" % sourceDb)
        for enzymeId, eD in dataD.items():
            self.assertEqual(eD["source_db"], sourceDb)
            self.assertTrue(eD["sequence"] and eD["name"], "%s %s is missing a sequence or name" % (sourceDb, enzymeId))
        self.assertTrue(any(eD["substrates"] for eD in dataD.values()), "no %s records have substrates" % sourceDb)
        self.assertTrue(peP.exportFasta(os.path.join(cachePath, "plastic-enzyme-sequences.fa")))

    def testPazySourceAccess(self):
        """The live PAZy source data is accessible and usable (fails if PAZy can't be fetched)."""
        self.__checkLiveSource("PAZy", "pazyUrl", "pazy-select-data.json")

    def testPlasticDBSourceAccess(self):
        """The live PlasticDB source data is accessible and usable (fails if PlasticDB can't be fetched)."""
        self.__checkLiveSource("PlasticDB", "plasticDBUrl", "plasticdb-select-data.json")


def suitePlasticEnzymeAnnotationProvider():
    suiteSelect = unittest.TestSuite()
    suiteSelect.addTest(PlasticEnzymeAnnotationProviderTests("testFetchAndStandardize"))
    suiteSelect.addTest(PlasticEnzymeAnnotationProviderTests("testFetchMinCount"))
    suiteSelect.addTest(PlasticEnzymeAnnotationProviderTests("testFetchFailure"))
    suiteSelect.addTest(PlasticEnzymeAnnotationProviderTests("testFetchFallbackToCache"))
    suiteSelect.addTest(PlasticEnzymeAnnotationProviderTests("testFetchFallbackToStash"))
    suiteSelect.addTest(PlasticEnzymeAnnotationProviderTests("testFetchUseCache"))
    suiteSelect.addTest(PlasticEnzymeAnnotationProviderTests("testTempLocalSourceFiles"))
    suiteSelect.addTest(PlasticEnzymeAnnotationProviderTests("testFetchStashKeepsCachedData"))
    suiteSelect.addTest(PlasticEnzymeAnnotationProviderTests("testExportFasta"))
    suiteSelect.addTest(PlasticEnzymeAnnotationProviderTests("testBuildAnnotationList"))
    suiteSelect.addTest(PlasticEnzymeAnnotationProviderTests("testReloadMalformedCache"))
    suiteSelect.addTest(PlasticEnzymeAnnotationProviderTests("testPazySourceAccess"))
    suiteSelect.addTest(PlasticEnzymeAnnotationProviderTests("testPlasticDBSourceAccess"))
    return suiteSelect


if __name__ == "__main__":
    mySuite = suitePlasticEnzymeAnnotationProvider()
    unittest.TextTestRunner(verbosity=2).run(mySuite)
