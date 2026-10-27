##
# File:    testPlasticEnzymeAnnotationProvider.py
# Date:    20-Oct-2026
#
# Updates:
#  26-Oct-2026 dwp Add testReloadMalformedCache
#  27-Oct-2026 dwp Replace the opt-in live fetch test with tests that always check that the PAZy and PlasticDB
#                  source data are accessible and usable (each fails if its source can't be fetched); start each test
#                  from an empty cache; add tests for the kept raw source files, the stash fallback (newest data
#                  kept), restore() using a valid local cache, and backup() adding missing raw source files; the live
#                  source tests also check that the downloaded raw source file is kept and not empty
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
from rcsb.utils.io.StashableBase import StashableBase
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
        # The raw source files are kept (for the stashed backup) next to the standardized data
        dirPath = os.path.join(self.__cachePath, "plastic-enzyme-annotations")
        self.assertEqual(sorted(os.listdir(dirPath)), ["PlasticDB.fasta", "pazy-select-data.json", "pazy_proteins.json", "plasticdb-select-data.json"])
        with open(os.path.join(dirPath, "PlasticDB.fasta"), encoding="utf-8") as ifh:
            self.assertEqual(ifh.read(), PLASTICDB_FASTA)
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
        """Source data with too few records isn't kept (neither the raw nor the standardized data)."""
        peP = self.__getTestProvider()
        ok = peP.fetchPlasticEnzymeData(minCount=100)
        self.assertFalse(ok)
        self.assertFalse(peP.testSourceCache(minCount=0))
        self.assertEqual(os.listdir(os.path.join(self.__cachePath, "plastic-enzyme-annotations")), [])

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

    def __writeStash(self, dirPath, tag, n=2):
        """Write a stand-in for a restored stash bundle (standardized data and raw files for both sources) to dirPath."""
        self.__mU.mkdir(dirPath)
        for fileName, sourceDb in [("pazy-select-data.json", "PAZy"), ("plasticdb-select-data.json", "PlasticDB")]:
            dataD = {"%s%d" % (tag, ii): {"source_db": sourceDb, "name": tag, "substrates": ["PET"], "sequence": SEQ_A} for ii in range(n)}
            self.__mU.doExport(os.path.join(dirPath, fileName), {"version": "1.0", "created": tag, "data": dataD}, fmt="json")
        for fileName in ["pazy_proteins.json", "PlasticDB.fasta"]:
            with open(os.path.join(dirPath, fileName), "w", encoding="utf-8") as ofh:
                ofh.write("%s %s\n" % (tag, fileName))

    def __writeLocalCache(self, n=101):
        """Write valid local standardized and annotation data (more than the default minCount of 100 records), without raw files."""
        dirPath = os.path.join(self.__cachePath, "plastic-enzyme-annotations")
        self.__writeStash(dirPath, "local", n=n)
        for fileName in ["pazy_proteins.json", "PlasticDB.fasta"]:
            os.remove(os.path.join(dirPath, fileName))
        aD = {"%dABC_1" % ii: [{"source_db": "PAZy", "annotation_id": "local0"}] for ii in range(n)}
        self.__mU.doExport(os.path.join(dirPath, "plastic-enzyme-annotation-data.json"), {"version": "1.0", "created": "local", "annotations": aD}, fmt="json")
        return dirPath

    def testFetchFallbackToStash(self):
        """If one source can't be fetched, the stashed backup is restored and the freshly fetched source is written over it."""
        dirPath = os.path.join(self.__cachePath, "plastic-enzyme-annotations")

        def fakeRestore(_, cfgOb, configName, remotePrefix=None, useStash=True, useGit=False):
            self.assertEqual((cfgOb, configName, remotePrefix, useStash, useGit), ("cfgOb", "configName", "A", True, False))
            self.__writeStash(dirPath, "stashed")
            return True

        peP = PlasticEnzymeAnnotationProvider(cachePath=self.__cachePath, useCache=False, pazyUrl=os.path.join(self.__workPath, "missing.json"), plasticDBUrl=self.__plasticDBPath)
        with mock.patch.object(StashableBase, "restore", autospec=True, side_effect=fakeRestore) as mRestore:
            self.assertFalse(peP.fetchPlasticEnzymeData(minCount=1, cfgOb="cfgOb", configName="configName", remotePrefix="A"))
        mRestore.assert_called_once()
        # PAZy comes from the stash (standardized and raw data)
        self.assertEqual(peP.getEnzyme("PAZy", "stashed0")["name"], "stashed")
        with open(os.path.join(dirPath, "pazy_proteins.json"), encoding="utf-8") as ifh:
            self.assertEqual(ifh.read(), "stashed pazy_proteins.json\n")
        # PlasticDB is the freshly fetched data, not the stashed copy (standardized and raw data)
        self.assertIsNotNone(peP.getEnzyme("PlasticDB", "00001"))
        self.assertIsNone(peP.getEnzyme("PlasticDB", "stashed0"))
        with open(os.path.join(dirPath, "PlasticDB.fasta"), encoding="utf-8") as ifh:
            self.assertEqual(ifh.read(), PLASTICDB_FASTA)
        self.assertTrue(peP.testSourceCache(minCount=1))
        self.assertFalse([fn for fn in os.listdir(dirPath) if fn.endswith(".download")])

    def testRestoreUsesValidLocalCache(self):
        """restore() uses a valid local cache (testCache()) instead of the stash, and restores from the stash otherwise."""
        dirPath = self.__writeLocalCache()
        peP = PlasticEnzymeAnnotationProvider(cachePath=self.__cachePath, useCache=False)
        with mock.patch.object(StashableBase, "restore", autospec=True, return_value=True) as mRestore:
            self.assertTrue(peP.restore("cfgOb", "configName", remotePrefix="A"))
        mRestore.assert_not_called()
        self.assertTrue(PlasticEnzymeAnnotationProvider(cachePath=self.__cachePath, useCache=True).testCache())
        # An incomplete local cache (no annotation data) is restored from the stash
        os.remove(os.path.join(dirPath, "plastic-enzyme-annotation-data.json"))
        with mock.patch.object(StashableBase, "restore", autospec=True, return_value=False) as mRestore:
            self.assertFalse(peP.restore("cfgOb", "configName", remotePrefix="A"))
        mRestore.assert_called_once()

    def testBackupAddsMissingRawFiles(self):
        """backup() adds missing raw source files from the stashed backup first, and skips the backup if they aren't there."""
        dirPath = self.__writeLocalCache()
        peP = PlasticEnzymeAnnotationProvider(cachePath=self.__cachePath, useCache=True)
        backedUpL = []

        def fakeRestore(tP, cfgOb, configName, remotePrefix=None, useStash=True, useGit=False):
            # Restore into the separate (temporary) cache directory of the provider instance it's called for
            _ = (cfgOb, configName, remotePrefix, useStash, useGit)
            self.__writeStash(os.path.join(tP._StashableBase__cachePath, "plastic-enzyme-annotations"), "stashed")  # pylint: disable=protected-access
            return True

        def fakeBackup(*_args, **_kwargs):
            backedUpL.append(sorted(os.listdir(dirPath)))
            return True

        with mock.patch.object(StashableBase, "restore", autospec=True, side_effect=fakeRestore), mock.patch.object(StashableBase, "backup", autospec=True, side_effect=fakeBackup):
            self.assertTrue(peP.backup("cfgOb", "configName", remotePrefix="A"))
        self.assertEqual(backedUpL, [sorted(["PlasticDB.fasta", "pazy_proteins.json", "pazy-select-data.json", "plasticdb-select-data.json", "plastic-enzyme-annotation-data.json"])])
        # Only the missing raw files are added; the local standardized data is unchanged
        self.assertEqual(peP.getEnzyme("PAZy", "local0")["name"], "local")
        self.assertEqual(self.__mU.doImport(os.path.join(dirPath, "pazy-select-data.json"), fmt="json")["created"], "local")
        self.assertFalse([fn for fn in os.listdir(self.__cachePath) if fn.startswith("plastic-enzyme-restore-")])
        #
        # Without the raw files in the stash, the backup is skipped
        for fileName in ["pazy_proteins.json", "PlasticDB.fasta"]:
            os.remove(os.path.join(dirPath, fileName))
        with mock.patch.object(StashableBase, "restore", autospec=True, return_value=False), mock.patch.object(StashableBase, "backup", autospec=True) as mBackup:
            self.assertFalse(peP.backup("cfgOb", "configName", remotePrefix="A"))
        mBackup.assert_not_called()

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

    def __checkLiveSource(self, sourceDb, urlKey, rawFileName, dataFileName, minCount=100):
        """Fetch one source from its default (live) URL, and the other from its local test file, and check that the live
        source's raw file is downloaded and kept, and returns enough usable records."""
        cachePath = os.path.join(self.__workPath, "CACHE-live-%s" % sourceDb)
        shutil.rmtree(cachePath, ignore_errors=True)
        kwargs = {"pazyUrl": self.__pazyPath, "plasticDBUrl": self.__plasticDBPath}
        kwargs.pop(urlKey)
        peP = PlasticEnzymeAnnotationProvider(cachePath=cachePath, useCache=False, **kwargs)
        ok = peP.fetchPlasticEnzymeData(minCount=1)
        self.assertTrue(ok, "%s source data could not be fetched (see the error logged above)" % sourceDb)
        rawPath = os.path.join(cachePath, "plastic-enzyme-annotations", rawFileName)
        self.assertTrue(os.path.isfile(rawPath) and os.path.getsize(rawPath) > 0, "%s raw source file %r is missing or empty" % (sourceDb, rawFileName))
        dataD = self.__mU.doImport(os.path.join(cachePath, "plastic-enzyme-annotations", dataFileName), fmt="json")["data"]
        self.assertGreater(len(dataD), minCount, "too few %s records" % sourceDb)
        for enzymeId, eD in dataD.items():
            self.assertEqual(eD["source_db"], sourceDb)
            self.assertTrue(eD["sequence"] and eD["name"], "%s %s is missing a sequence or name" % (sourceDb, enzymeId))
        self.assertTrue(any(eD["substrates"] for eD in dataD.values()), "no %s records have substrates" % sourceDb)
        self.assertTrue(peP.exportFasta(os.path.join(cachePath, "plastic-enzyme-sequences.fa")))

    def testPazySourceAccess(self):
        """The live PAZy source data is accessible and usable (fails if PAZy can't be fetched)."""
        self.__checkLiveSource("PAZy", "pazyUrl", "pazy_proteins.json", "pazy-select-data.json")

    def testPlasticDBSourceAccess(self):
        """The live PlasticDB source data is accessible and usable (fails if PlasticDB can't be fetched)."""
        self.__checkLiveSource("PlasticDB", "plasticDBUrl", "PlasticDB.fasta", "plasticdb-select-data.json")


def suitePlasticEnzymeAnnotationProvider():
    suiteSelect = unittest.TestSuite()
    suiteSelect.addTest(PlasticEnzymeAnnotationProviderTests("testFetchAndStandardize"))
    suiteSelect.addTest(PlasticEnzymeAnnotationProviderTests("testFetchMinCount"))
    suiteSelect.addTest(PlasticEnzymeAnnotationProviderTests("testFetchFailure"))
    suiteSelect.addTest(PlasticEnzymeAnnotationProviderTests("testFetchFallbackToStash"))
    suiteSelect.addTest(PlasticEnzymeAnnotationProviderTests("testRestoreUsesValidLocalCache"))
    suiteSelect.addTest(PlasticEnzymeAnnotationProviderTests("testBackupAddsMissingRawFiles"))
    suiteSelect.addTest(PlasticEnzymeAnnotationProviderTests("testExportFasta"))
    suiteSelect.addTest(PlasticEnzymeAnnotationProviderTests("testBuildAnnotationList"))
    suiteSelect.addTest(PlasticEnzymeAnnotationProviderTests("testReloadMalformedCache"))
    suiteSelect.addTest(PlasticEnzymeAnnotationProviderTests("testPazySourceAccess"))
    suiteSelect.addTest(PlasticEnzymeAnnotationProviderTests("testPlasticDBSourceAccess"))
    return suiteSelect


if __name__ == "__main__":
    mySuite = suitePlasticEnzymeAnnotationProvider()
    unittest.TextTestRunner(verbosity=2).run(mySuite)
