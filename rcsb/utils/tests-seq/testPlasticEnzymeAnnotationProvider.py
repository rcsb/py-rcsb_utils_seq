##
# File:    testPlasticEnzymeAnnotationProvider.py
# Date:    20-Oct-2026
#
# Updates:
#
##
"""
Tests for the plastic-degrading enzyme (PAZy and PlasticDB) annotation provider.
"""

import json
import logging
import os
import unittest

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

    @unittest.skipUnless(os.environ.get("PLASTIC_ENZYME_LIVE_TEST"), "Set PLASTIC_ENZYME_LIVE_TEST to fetch the live PAZy and PlasticDB data")
    def testFetchLive(self):
        """Fetch the live PAZy and PlasticDB data."""
        peP = PlasticEnzymeAnnotationProvider(cachePath=os.path.join(self.__workPath, "CACHE-live"), useCache=False)
        self.assertTrue(peP.fetchPlasticEnzymeData())
        self.assertTrue(peP.testSourceCache())
        self.assertTrue(peP.exportFasta())


def suitePlasticEnzymeAnnotationProvider():
    suiteSelect = unittest.TestSuite()
    suiteSelect.addTest(PlasticEnzymeAnnotationProviderTests("testFetchAndStandardize"))
    suiteSelect.addTest(PlasticEnzymeAnnotationProviderTests("testFetchMinCount"))
    suiteSelect.addTest(PlasticEnzymeAnnotationProviderTests("testFetchFailure"))
    suiteSelect.addTest(PlasticEnzymeAnnotationProviderTests("testExportFasta"))
    suiteSelect.addTest(PlasticEnzymeAnnotationProviderTests("testBuildAnnotationList"))
    suiteSelect.addTest(PlasticEnzymeAnnotationProviderTests("testFetchLive"))
    return suiteSelect


if __name__ == "__main__":
    mySuite = suitePlasticEnzymeAnnotationProvider()
    unittest.TextTestRunner(verbosity=2).run(mySuite)
