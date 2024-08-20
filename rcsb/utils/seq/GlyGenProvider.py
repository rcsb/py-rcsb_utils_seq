##
#  File:  GlyGenProvider.py
#  Date:  26-May-2021 jdw
#
#  Updates:
#  15-May-2022 dwp Update resource URL to new location
#  14-Nov-2023 dwp Update GlyGen data version and add version information to cache file;
#                  Add functionality for fetching data via SPARQL
#  20-Aug-2024 dwp Adjust fetch method following certificate changes
##
"""
  Fetch glycans and glycoproteins available in the GlyGen.org resource.

"""

import logging
import os.path
from SPARQLWrapper import SPARQLWrapper, JSON
from rcsb.utils.io.FileUtil import FileUtil
from rcsb.utils.io.MarshalUtil import MarshalUtil
from rcsb.utils.io.StashableBase import StashableBase
from rcsb.utils.io.UrlRequestUtil import UrlRequestUtil

logger = logging.getLogger(__name__)


class GlyGenProvider(StashableBase):
    """Fetch glycans and glycoproteins available in the GlyGen.org resource.

    GlyGen glycan link template -
          https://glygen.org/glycan/G28882EF

    Glycoprotein link template -
          https://www.glygen.org/protein/Q658T7
    """

    def __init__(self, **kwargs):
        #
        dirName = "glygen"
        cachePath = kwargs.get("cachePath", ".")
        self.__dirPath = os.path.join(cachePath, dirName)
        super(GlyGenProvider, self).__init__(cachePath, [dirName])
        useCache = kwargs.get("useCache", True)
        #
        baseUrl = kwargs.get("glygenBasetUrl", "https://data.glygen.org/ln2data/releases/data/v-2.2.1/reviewed/")
        # baseSparqlUrl = kwargs.get("glygenBaseSparqlUrl", "http://sparql.glygen.org:8880/sparql")
        fallbackUrl = kwargs.get("glygenFallbackUrl", "https://raw.githubusercontent.com/rcsb/py-rcsb_exdb_assets/master/fall_back/glygen/")
        #
        self.__mU = MarshalUtil(workPath=self.__dirPath)
        self.__glycanD, self.__version = self.__reloadGlycans(baseUrl, fallbackUrl, self.__dirPath, useCache=useCache)
        self.__glycoproteinD = self.__reloadGlycoproteins(baseUrl, fallbackUrl, self.__dirPath, useCache=useCache)
        # self.__glycoproteinD = self.__reloadGlycoproteinsSparql(baseSparqlUrl, fallbackUrl, self.__dirPath, useCache=useCache)

    def testCache(self, minGlycanCount=20000, minGlycoproteinCount=64000):
        #
        logger.info("GlyGen glycan list (%d) glycoprotein list (%d)", len(self.__glycanD), len(self.__glycoproteinD))
        if self.__glycanD and len(self.__glycanD) > minGlycanCount and self.__glycoproteinD and len(self.__glycoproteinD) > minGlycoproteinCount:
            return True
        return False

    def hasGlycan(self, glyTouCanId):
        try:
            return glyTouCanId in self.__glycanD
        except Exception:
            return False

    def hasGlycoprotein(self, uniProtId):
        try:
            return uniProtId in self.__glycoproteinD
        except Exception:
            return False

    def getVersion(self):
        return self.__version

    def getGlycans(self):
        return self.__glycanD

    def getGlycoproteins(self):
        return self.__glycoproteinD

    def __reloadGlycans(self, baseUrl, fallbackUrl, dirPath, useCache=True):
        gD = {}
        version = "1.0"
        logger.debug("Using dirPath %r", dirPath)
        self.__mU.mkdir(dirPath)
        #
        myDataPath = os.path.join(dirPath, "glygen-glycan-list.json")
        if useCache and self.__mU.exists(myDataPath):
            fD = self.__mU.doImport(myDataPath, fmt="json")
            gD = fD["data"]
            version = fD["version"]
            logger.debug("GlyGen glycan data length %d", len(gD))
        elif not useCache:
            logger.debug("Fetch GlyGen glycan data from primary data source %s", baseUrl)
            endPoint = os.path.join(baseUrl, "glycan_masterlist.csv")
            #
            logger.info("Fetch GlyGen glycan data from primary data source %s", endPoint)
            rawPath = os.path.join(dirPath, "glycan_masterlist.csv")
            fU = FileUtil()
            uR = UrlRequestUtil()
            ret, retCode = uR.get(baseUrl, "glycan_masterlist.csv", {})
            ok = retCode == 200
            logger.debug("Fetch GlyGen glycan data status %r", ok)
            if ok:
                with open(rawPath, "w", encoding="utf-8") as f:
                    f.write(ret)
            #
            versionEndPoint = os.path.join(baseUrl, "release-notes.txt")
            try:
                uR = UrlRequestUtil()
                ret, retCode = uR.get(baseUrl, "release-notes.txt", {})
                version = ret.split(" ")[0].split("v-")[-1]
            except Exception as e:
                logger.exception("Failing for %r with %s", versionEndPoint, str(e))
            #
            if not ok:
                endPoint = os.path.join(fallbackUrl, "glycan_masterlist.csv")
                ok = fU.get(endPoint, rawPath)
                version = "1.0"
                logger.info("Fetch fallback GlyGen glycan data status %r", ok)
            #
            if ok:
                gD = self.__parseGlycanList(rawPath)
                fD = {"data": gD, "version": version}
                ok = self.__mU.doExport(myDataPath, fD, fmt="json")
                logger.info("Exported GlyGen glycan list (%d) version (%r) (%r) %s", len(gD), version, ok, myDataPath)
            #
        return gD, version

    def __parseGlycanList(self, filePath):
        gD = {}
        row = None
        try:
            rowL = self.__mU.doImport(filePath, fmt="csv", rowFormat="list")
            logger.debug("Glycan list length (%d)", len(rowL))
            logger.debug("Row 0 %r", rowL[0])
            for row in rowL[1:]:
                gD[row[0]] = row[1]
        except Exception as e:
            logger.exception("Failing for %r (%r) with %s", filePath, row, str(e))
        return gD

    def __reloadGlycoproteins(self, baseUrl, fallbackUrl, dirPath, useCache=True):
        gD = {}
        version = "1.0"
        logger.debug("Using dirPath %r", dirPath)
        self.__mU.mkdir(dirPath)
        #
        myDataPath = os.path.join(dirPath, "glygen-glycoprotein-list.json")
        if useCache and self.__mU.exists(myDataPath):
            fD = self.__mU.doImport(myDataPath, fmt="json")
            gD = fD["data"]
            version = fD["version"]
            logger.debug("GlyGen glycoprotein data length %d", len(gD))
        else:
            #
            versionEndPoint = os.path.join(baseUrl, "release-notes.txt")
            try:
                uR = UrlRequestUtil()
                ret, retCode = uR.get(baseUrl, "release-notes.txt", {})
                version = ret.split(" ")[0].split("v-")[-1]
            except Exception as e:
                logger.exception("Failing for %r with %s", versionEndPoint, str(e))
            #
            for fn in [
                "sarscov1_protein_masterlist.csv",
                "sarscov2_protein_masterlist.csv",
                "hcv1b_protein_masterlist.csv",
                "hcv1a_protein_masterlist.csv",
                "human_protein_masterlist.csv",
                "mouse_protein_masterlist.csv",
                "rat_protein_masterlist.csv",
                "fruitfly_protein_masterlist.csv",
                "yeast_protein_masterlist.csv",
            ]:
                logger.debug("Fetch GlyGen glycoprotein data from primary data source %s", baseUrl)
                endPoint = os.path.join(baseUrl, fn)
                #
                logger.debug("Fetch GlyGen glycoprotein data from primary data source %s", endPoint)
                rawPath = os.path.join(dirPath, fn)
                fU = FileUtil()
                uR = UrlRequestUtil()
                ret, retCode = uR.get(baseUrl, fn, {})
                ok = retCode == 200
                logger.info("Fetch GlyGen glycoprotein data status %r - %r", ok, endPoint)
                if ok:
                    with open(rawPath, "w", encoding="utf-8") as f:
                        f.write(ret)
                else:
                    # Fetch from fallback
                    endPoint = os.path.join(fallbackUrl, fn)
                    ok = fU.get(endPoint, rawPath)
                    logger.info("Fetch fallback GlyGen data status %r - %r", ok, endPoint)
                    version = "1.0"
                #
                if ok:
                    tD = self.__parseGlycoproteinList(rawPath)
                    gD.update(tD)
            #
            fD = {"data": gD, "version": version}
            ok = self.__mU.doExport(myDataPath, fD, fmt="json")
            logger.info("Exported GlyGen glycoprotein list (%d) version (%r) (%r) %s", len(gD), version, ok, myDataPath)
        #
        return gD

    def __parseGlycoproteinList(self, filePath):
        gD = {}
        try:
            rowL = self.__mU.doImport(filePath, fmt="csv", rowFormat="list")
            for row in rowL[1:]:
                ff = row[0].split("-")
                gD[ff[0]] = ff[1]
        except Exception as e:
            logger.exception("Failing for %r with %s", filePath, str(e))
        return gD

    def __reloadGlycoproteinsSparql(self, baseSparqlUrl, dirPath, useCache=True):
        gD = {}
        logger.debug("Using dirPath %r", dirPath)
        self.__mU.mkdir(dirPath)
        #
        myDataPath = os.path.join(dirPath, "glygen-glycoprotein-list.json")
        if useCache and self.__mU.exists(myDataPath):
            gD = self.__mU.doImport(myDataPath, fmt="json")
            logger.info("GlyGen glycoprotein data length %d", len(gD))
        else:
            for organism, taxId in {
                "sarscov1": "694009",
                "sarscov2": "2697049",
                "hcv1b": "11116",
                "hcv1a": "63746",
                "human": "9606",
                "mouse": "10090",
                "rat": "10116",
                "fruitfly": "7227",
                "yeast": "4932",
            }.items():
                logger.info("Fetch GlyGen glycoprotein data for organism %s taxId %s from SPARQL source %s", organism, taxId, baseSparqlUrl)
                resultL, retL = [], []
                offset = 0
                retL = self.__fetchGlycoproteinListSparql(baseSparqlUrl, taxId, offset)
                while len(retL) > 0:
                    resultL += retL
                    offset += len(retL)
                    retL = self.__fetchGlycoproteinListSparql(baseSparqlUrl, taxId, offset)
                logger.info("GlyGen glycoprotein data length (%d) for organism %s taxId %s", len(resultL), organism, taxId)

                if len(resultL) > 0:
                    tD = {}
                    for r in resultL:
                        ff = r.split("-")
                        tD[ff[0]] = ff[1]
                    gD.update(tD)

            # if len(gD) < 100:
            #     endPoint = os.path.join(fallbackUrl, fn)
            #     ok = fU.get(endPoint, rawPath)
            #     logger.info("Fetch fallback GlyGen data status %r", ok)

            ok = self.__mU.doExport(myDataPath, gD, fmt="json")
            logger.info("Exported GlyGen glycoprotein list (%d) (%r) %s", len(gD), ok, myDataPath)
        return gD

    def __fetchGlycoproteinListSparql(self, baseSparqlUrl, taxId, offset):
        retL = []
        try:
            # Set the SPARQL endpoint
            sparql = SPARQLWrapper(baseSparqlUrl)

            # Define the SPARQL query
            query = f"""
            PREFIX glycan: <http://purl.jp/bio/12/glyco/glycan#>
            PREFIX up: <http://purl.uniprot.org/core/>
            PREFIX gly:<https://sparql.glygen.org/ontology/>
            PREFIX skos:<http://www.w3.org/2004/02/skos/core#>
            SELECT DISTINCT ?isoform_uri
            WHERE {{
                ?glycoprotein_uri up:sequence ?isoform_uri .
                ?protein_uri up:sequence ?isoform_uri .
                ?protein_uri up:organism <http://purl.uniprot.org/taxonomy/{taxId}> .
                ?isoform_uri gly:canonical "true"^^<http://www.w3.org/2001/XMLSchema#boolean> .
            }}
            LIMIT 10000 OFFSET {offset}
            """

            # Set the query and the return format
            sparql.setQuery(query)
            sparql.setReturnFormat(JSON)

            # Execute the query
            results = sparql.query().convert()

            # Process the results
            for result in results["results"]["bindings"]:
                isoformUri = result["isoform_uri"]["value"]
                retL.append(isoformUri.split("/")[-1])

        except Exception as e:
            logger.exception("Failing for taxId %s offset with %s with %s", taxId, offset, str(e))

        return retL
