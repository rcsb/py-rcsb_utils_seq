##
# File:    UniProtUtils.py
# Date:    14-Mar-2019
#
# JDW - Adapted from SBKB sources and bits from RCSB codes -
#
# Updates:
#  6-Dec-2019 jdw Add method rebuildMatchResultIndex() to refresh match index from existing reference store
# 25-Jul-2022 dwp Adjust UniProt API calls to temporarily use the legacy baseUrl, legacy.uniprot.org
#  4-Oct-2022 dwp Update code to use new UniProt API (much of the code from: https://www.uniprot.org/help/id_mapping)
# 11-Oct-2022 dwp Only use secondary service site for bulk XML fetchList requests
# 26-Nov-2022 dwp Fix primary service fetching method (use accessions endpoint, not id_mapping)
# 28-Nov-2022 dwp Add functionality to search primary service using secondary accession IDs (in case of obsoleted IDs)
# 23-Feb-2023 dwp Fix primary service fetching method to handle failures due to one or more invalid UniProt IDs
#  9-Mar-2023 dwp Lower default maxChunkSize to 10 (UniProt API having trouble streaming XML responses)
#
##

import collections
import json
import logging
import os
import time
import re
import requests

from rcsb.utils.io.FastaUtil import FastaUtil
from rcsb.utils.io.UrlRequestUtil import UrlRequestUtil
from rcsb.utils.seq.UniProtReader import UniProtReader

try:
    from itertools import zip_longest
except ImportError:
    # Python 2
    from itertools import izip_longest as zip_longest

logger = logging.getLogger(__name__)

FeatureLabel = collections.namedtuple("FeatureLabel", "type description id evidence reference original variation")


class UniProtUtils(object):
    """
    Manage fetch queries for UniProt entries and related annotations.

    XML entry data is parsed into a feature dictionary.
    """

    def __init__(self, **kwargs):
        self.__saveText = kwargs.get("saveText", False)
        self.__urlPrimary = kwargs.get("urlPrimary", "https://rest.uniprot.org")
        self.__urlSecondary = kwargs.get("urlSecondary", "https://www.ebi.ac.uk")
        self.__dataList = []
        self.__unpFeatureD = {
            kk: kk
            for kk in [
                "ACTIVE_SITE",
                "BINDING_SITE",
                "CALCIUM_BINDING_REGION",
                "CHAIN",
                "COILED_COIL_REGION",
                "COMPOSITIONALLY_BIASED_REGION",
                "CROSS_LINK",
                "DISULFIDE_BOND",
                "DNA_BINDING_REGION",
                "DOMAIN",
                "GLYCOSYLATION_SITE",
                "HELIX",
                "INITIATOR_METHIONINE",
                "LIPID_MOIETY_BINDING_REGION",
                "METAL_ION_BINDING_SITE",
                "MODIFIED_RESIDUE",
                "MUTAGENESIS_SITE",
                "NON_CONSECUTIVE_RESIDUES",
                "NON_TERMINAL_RESIDUE",
                "NUCLEOTIDE_PHOSPHATE_BINDING_REGION",
                "PEPTIDE",
                "PROPEPTIDE",
                "REGION_OF_INTEREST",
                "REPEAT",
                "NON_STANDARD_AMINO_ACID",
                "SEQUENCE_CONFLICT",
                "SEQUENCE_VARIANT",
                "SHORT_SEQUENCE_MOTIF",
                "SIGNAL_PEPTIDE",
                "SITE",
                "SPLICE_VARIANT",
                "STRAND",
                "TOPOLOGICAL_DOMAIN",
                "TRANSIT_PEPTIDE",
                "TRANSMEMBRANE_REGION",
                "TURN",
                "UNSURE_RESIDUE",
                "ZINC_FINGER_REGION",
                "INTRAMEMBRANE_REGION",
            ]
        }
        #

    def fetchList(self, idList, maxChunkSize=10, usePrimary=True, retryAltApi=True):
        """Execute a fetch query for the input id list. The input list is filtered
           for sequence variants (e.g. ids with appended '-#').

           The input list is divided into manageable chunks, each chunk is fetched,
           and result concatenated.

        Args:
            idList (list): List of UniProt identifiers
            maxChunkSize (int, optional): id chunk size. Defaults to 100.
            usePrimary (bool, optional): Use the primary UniProt webservice. Defaults to True.

        Returns:
            dict: {unpId: {'key':val, ... }} dictionary of UniProt reference data
            dict: {unpId: {match details}, ...} } dictionary of match details
        """
        try:
            if self.__saveText:
                self.__dataList = []
            referenceD = {}
            matchD = {}
            #
            searchIdList, variantD = self.__processIdList(idList)

            logger.info("Running fetchList for idList length %d with maxChunkSize %d usePrimary %r retryAltApi %r", len(idList), maxChunkSize, usePrimary, retryAltApi)
            logger.debug("input id list %s", idList)
            logger.debug("search   list %s", searchIdList)
            logger.debug("variants      %s", variantD.items())

            if searchIdList is None or not searchIdList:
                return referenceD, matchD
            #
            searchIdList = list(set(searchIdList))
            subLists = self.__makeSubLists(maxChunkSize, searchIdList)
            # numLists = len(searchIdList) / maxChunkSize + 1
            # numLists = len(subLists)

            for ii, subList in enumerate(subLists):
                missingIds, demergedIdList = [], []
                logger.debug("Fetching subList %r", subList)
                logger.debug("Starting fetching for sublist %d", ii + 1)
                #
                ok, xmlText, missingIds, demergedIdList = self.__doRequest(subList, usePrimary=usePrimary, retryAltApi=retryAltApi)
                if not ok and xmlText is not None:
                    logger.error("Failing request with status %r, len(xmlText) %r, xmlText[-100:] %r", ok, len(xmlText), xmlText[-100:])
                logger.debug("Status %r", ok)
                #
                # Filter possible simple text error messages from the failed queries.
                #
                if (xmlText is not None) and not xmlText.startswith("ERROR"):
                    tD = self.__parseText(xmlText, variantD)
                    if tD:
                        referenceD.update(tD)
                    else:
                        logger.error("Status %r. Bad xml text. First 100 characters: %r; Last 100 characters: %r ", ok, xmlText[:100], xmlText[-100:])
                        raise ValueError("Bad xml text")
                    if self.__saveText:
                        self.__dataList.append(xmlText)
                else:
                    logger.info("Fetch %r status %r text %r missingIds %r demergedIdList %r", subList, ok, xmlText, missingIds, demergedIdList)

                # Re-attempt fetch for missingIds (by searching for corresponding demergedIdList -- i.e., secondary accession ids)
                if len(demergedIdList) > 0:
                    logger.info("Re-fetching missingIds %r using demergedIdList %r", missingIds, demergedIdList)
                    ok, xmlText, missingIds2, demergedIdList2 = self.__doRequest(demergedIdList, usePrimary=usePrimary, retryAltApi=retryAltApi)
                    if not ok and xmlText is not None:
                        logger.error("Failing request with status %r, len(xmlText) %r, xmlText[-100:] %r", ok, len(xmlText), xmlText[-100:])
                    logger.debug("Status %r", ok)
                    #
                    # Filter possible simple text error messages from the failed queries.
                    #
                    if (xmlText is not None) and not xmlText.startswith("ERROR"):
                        tD = self.__parseText(xmlText, variantD)
                        if tD:
                            referenceD.update(tD)
                        else:
                            logger.error("Status %r. Bad xml text. First 100 characters: %r; Last 100 characters: %r ", ok, xmlText[:100], xmlText[-100:])
                            raise ValueError("Bad xml text")
                        if self.__saveText:
                            self.__dataList.append(xmlText)
                    else:
                        logger.info("Fetch %r status %r text %r missingIds %r demergedIdList %r", subList, ok, xmlText, missingIds2, demergedIdList2)

            matchD = self.rebuildMatchResultIndex(idList, referenceD)

        except Exception as e:
            logger.exception("Failing with %s", str(e))

        return referenceD, matchD

    def rebuildMatchResultIndex(self, idList, referenceD):
        """Rebuild the search result index for the input id list.
        The input list is filtered for variants (e.g. ids with appended '-#').

        Divide the input list into manageable chunks, fetch each chunk,
        and concatenate the result.

        Return dict: dictionary index of search match results

        """
        try:
            matchD = {}
            #
            searchIdList, variantD = self.__processIdList(idList)

            logger.debug("input id list %s", idList)
            logger.debug("search   list %s", searchIdList)
            logger.debug("variants      %s", variantD.items())

            if searchIdList is None or not searchIdList:
                return matchD
            #
            matchD = {inpId: {"searchId": searchIdList[ii]} for ii, inpId in enumerate(idList)}
            searchIdList = list(set(searchIdList))
            #
            # Create a match dictionary for the input id list -
            #
            for _, sD in matchD.items():
                sId = sD["searchId"]
                if "matched" in sD:
                    continue
                if sId in referenceD:
                    taxId = referenceD[sId]["taxonomy_id"] if "taxonomy_id" in referenceD[sId] else None
                    sD.setdefault("matchedIds", {}).update({sId: {"taxId": taxId}})
                    sD["matched"] = "primary"
                else:
                    for _, rD in referenceD.items():
                        if sId in rD["accessions"]:
                            pId = rD["db_accession"]
                            taxId = referenceD[pId]["taxonomy_id"] if "taxonomy_id" in referenceD[pId] else None
                            sD.setdefault("matchedIds", {}).update({pId: {"taxId": taxId}})
                            sD["matched"] = "secondary"
                    #
                    if "matched" not in sD:
                        sD["matched"] = "none"

        except Exception as e:
            logger.exception("Failing with %s", str(e))

        return matchD

    def reformat(self, refD, formatType="exchange"):
        defAssertion = "ECO:0000323"
        rObj = {}
        if formatType != "exchange":
            return rObj
        #
        for uId, uD in refD.items():
            evD = uD["evidence"] if "evidence" in uD else {}
            rD = {}
            dbName = uD["db_name"]
            dbVersion = uD["version"]
            rD["rcsb_id"] = uId
            rD["rcsb_uniprot_container_identifiers"] = {"uniprot_id": uId}
            if "accessions" in uD:
                rD["rcsb_uniprot_accession"] = uD["accessions"]
            if "db_code" in uD:
                rD.setdefault("rcsb_uniprot_entry_name", []).append(uD["db_code"])
            if "keywords" in uD:
                rD["rcsb_uniprot_keyword"] = [{"id": tD["id"], "value": tD["keyword"]} for tD in uD["keywords"]]
            #
            if "sequence" in uD:
                rD["rcsb_uniprot_protein"] = {"sequence": uD["sequence"]}
            if "names" in uD:
                for nD in uD["names"]:
                    if "nameType" in nD and nD["nameType"] in ["recommendedName"]:
                        rD["rcsb_uniprot_protein"]["name"] = {"value": nD["name"], "provenance_code": defAssertion}
            #
            if "gene" in uD:
                tL = []
                for nD in uD["gene"]:
                    # 	"enum" : [  "PRIMARY", "SYNONYM",  "ORDERED_LOCUS", "ORF"]
                    tt = nD["type"].upper().replace(" ", "_")
                    tL.append({"type": tt, "value": nD["name"]})
                rD["rcsb_uniprot_protein"].setdefault("gene", []).append({"name": tL})
            #
            if "comments" in uD:
                for cD in uD["comments"]:
                    if "type" in cD and cD["type"] == "function" and "text" in cD:
                        evKy = cD["evidence"] if "evidence" in cD and cD["evidence"] else None
                        evC = evD[evKy] if evKy in evD else defAssertion
                        txt = cD["text"]
                        rD["rcsb_uniprot_protein"]["function"] = {"details": txt, "provenance_code": evC}
            if "source_scientific" in uD and "taxonomy_id" in uD:
                evKy = uD["taxonomy_evc"] if "taxonomy_evc" in uD and uD["taxonomy_evc"] else None
                evC = evD[evKy] if evKy in evD else defAssertion
                rD["rcsb_uniprot_protein"]["source_organism"] = {"scientific_name": uD["source_scientific"], "taxonomy_id": uD["taxonomy_id"], "provenance_code": evC}
            #
            pfamL = set()
            ensL = set()
            goL = set()
            ecD = {}
            if "dbReferences" in uD:
                for tD in uD["dbReferences"]:
                    evKy = tD["evidence"] if "evidence" in tD and tD["evidence"] else None
                    evCode = evD[evKy] if evKy in evD else defAssertion
                    rsName = tD["resource"] if "resource" in tD and tD["resource"] else None
                    idCode = tD["id_code"] if "id_code" in tD and tD["id_code"] else None
                    if rsName == "EC":
                        ecD[idCode] = evCode
                    elif rsName == "Pfam":
                        pfamL.add(idCode)
                    elif rsName == "GO":
                        goL.add(idCode)
                    elif rsName.upper().startswith("ENSEMB"):
                        ensL.add(idCode)
            if pfamL:
                rD["rcsb_uniprot_container_identifiers"]["pfam_ids"] = sorted(pfamL)
            if goL:
                rD["rcsb_uniprot_container_identifiers"]["go_ids"] = sorted(goL)
            if ensL:
                rD["rcsb_uniprot_container_identifiers"]["ensembl_ids"] = sorted(ensL)
            if ecD:
                rD["rcsb_uniprot_protein"]["ec"] = [{"number": ecId, "provenance_code": pC} for ecId, pC in ecD.items()]
            # -------
            # create index of features =
            # FeatureLabel = collections.namedtuple("FeatureLabel", "type description id evidence reference orginal variation")
            #
            fIndx = {}
            if "features" in uD:
                fIndxD = {}
                for fD in uD["features"]:
                    fType = fD["type"].upper().replace(" ", "_")
                    fDes = fD["description"] if "description" in fD else None
                    fId = fD["feature_id"] if "feature_id" in fD else None
                    fRef = fD["reference"] if "reference" in fD else None
                    fEv = fD["evidence"] if "evidence" in fD else None
                    fOrg = fD["original"] if "original" in fD else None
                    fVar = fD["variation"] if "variation" in fD else None
                    fLabel = FeatureLabel(fType, fDes, fId, fEv, fRef, fOrg, fVar)
                    fIndxD.setdefault(fLabel, []).append(fD)
                #
                for fIdx, fDL in fIndxD.items():
                    if fIdx.type not in self.__unpFeatureD:
                        continue
                    tD = {}
                    tD["type"] = fIdx.type
                    tD["assignment_version"] = dbName + "_" + dbVersion
                    if fIdx.id:
                        tD["feature_id"] = fIdx.id
                    if fIdx.evidence:
                        try:
                            evL = [evD[t] for t in fIdx.evidence.split()]
                            tD["provenance_code"] = ",".join(evL)
                        except Exception as e:
                            logger.exception("Failing with %s", str(e))
                    #
                    tD["reference_scheme"] = "UniProt"
                    dS = ""
                    if fIdx.description:
                        dS = fIdx.description
                    if fIdx.original and fIdx.variation:
                        dS += " (" + fIdx.original + " -> " + fIdx.variation + ")"
                    if dS:
                        tD["description"] = dS
                    #
                    for fD in fDL:
                        if "begin" in fD and "end" in fD:
                            tD.setdefault("feature_ranges", []).append({"beg_seq_id": int(fD["begin"]), "end_seq_id": int(fD["end"])})
                        elif "position" in fD:
                            if fIdx.original:
                                tD.setdefault("feature_positions", []).append({"comp_id": fIdx.original, "seq_id": int(fD["position"])})
                            else:
                                tD.setdefault("feature_positions", []).append({"seq_id": int(fD["position"])})
                    #
                    serialD = json.dumps(tD, sort_keys=True)
                    if serialD in fIndx:
                        continue
                    fIndx[serialD] = True
                    rD.setdefault("rcsb_uniprot_feature", []).append(tD)
            #
            rObj[uId] = rD
        #
        return rObj

    def writeUnpXml(self, filePath):
        with open(filePath, "w", encoding="utf-8") as ofh:
            for data in self.__dataList:
                ofh.write(data)

    def __parseText(self, xmlText, variantD):
        """
        Parse the accumulated xml text for each chunk and store the parsed data in
        the internal result dictionary.

        variant id codes are regisistered with the parser so that the returned dictionary
        contains keys corresponding to the original input id list.
        """
        ur = UniProtReader()
        try:
            logger.debug("variantD %r", variantD)
            for vId, aId in variantD.items():
                ur.addVariant(aId, vId)
            retD = ur.readString(xmlText)
            logger.debug("retD: %r", retD)
            return retD
        except Exception as e:
            logger.exception("Failing with %s", str(e))
        return {}

    def __doRequest(self, idList, retryAltApi=True, usePrimary=True):
        ok = False
        missingIds, demergedIdList = [], []
        primaryCounter = 0
        #
        if usePrimary:
            while not ok and primaryCounter < 5:
                primaryCounter += 1
                logger.info("Performing primary request for idList length %d (attempt %d)", len(idList), primaryCounter)
                ret, retCode, missingIds, demergedIdList = self.__doRequestPrimary(idList)
                if retCode in [200] and ret is not None:
                    ok = len(ret) > 0 and "ERROR" not in ret[0:100].upper() and "ERROR" not in ret[-100:].upper()
                    logger.info("  success %r", ok)
                    if not ok:
                        logger.info("Beginning of response: %r", ret[0:100])
                        logger.info("      End of response: %r", ret[-100:])
        #
        if retryAltApi and not ok:
            logger.info("Retrying using secondary service site")
            ret, retCode = self.__doRequestSecondary(idList)
            ok = retCode in [200] and ret and len(ret) > 0
        #
        return ok, ret, missingIds, demergedIdList

    def __doRequestPrimary(self, idList):
        """ """
        missingIds, demergedIdList = [], []
        baseUrl = self.__urlPrimary
        ureq = UrlRequestUtil()
        idListStr = "%2C%20".join(idList)
        endPoint = "uniprotkb/accessions?accessions=" + idListStr
        hD = {"Accept": "application/xml"}
        logger.debug("Performing request with baseUrl, endpoint, and hD:  %r  %r  %r", baseUrl, endPoint, hD)
        ret, respCode = ureq.getUnWrapped(baseUrl, endPoint, paramD=None, headers=hD, overwriteUserAgent=False)
        if respCode != 200:
            logger.error("Primary request failed with retCode %r", respCode)
            res = requests.get(os.path.join(baseUrl, endPoint), timeout=600)
            if res.status_code == 400:
                if "It should be a valid UniProtKB accession" in res.json()["messages"][0]:
                    for msg in res.json()["messages"]:
                        badCode = msg.split("has invalid format")[0].split("Accession")[1].split("'")[1].strip()
                        missingIds.append(badCode)
                        idList.pop(idList.index(badCode))
                    idListStr = "%2C%20".join(idList)
                    endPoint = "uniprotkb/accessions?accessions=" + idListStr
                    logger.info("Retrying without badCodes %r", missingIds)
                    res = requests.get(os.path.join(baseUrl, endPoint), timeout=600)
            if res.status_code == 200:
                ret, respCode = ureq.getUnWrapped(baseUrl, endPoint, paramD=None, headers=hD, overwriteUserAgent=False)
        #
        for accId in idList:
            if "<accession>" + accId + "</accession>" not in ret:
                missingIds.append(accId)
                try:
                    endpoint2 = os.path.join("uniprotkb", accId)
                    ret2, respCode2 = ureq.getUnWrapped(baseUrl, endpoint2, paramD=None, headers={"Accept": "application/json"}, overwriteUserAgent=False, returnContentType="JSON")
                    if respCode2 == 200:
                        if "inactiveReason" in ret2:
                            demergedIds = ret2["inactiveReason"].get("mergeDemergeTo", [])
                            if len(demergedIds) > 0:
                                demergedIdList.extend(ret2["inactiveReason"]["mergeDemergeTo"])
                except Exception as e:
                    logger.exception("Failing to collect demergedIdList with message %s", str(e))

        return ret, respCode, missingIds, demergedIdList

    def __doRequestSecondary(self, idList):
        baseUrl = self.__urlSecondary
        endPoint = "proteins/api/proteins"
        #
        hL = [("Accept", "application/xml")]
        pD = {}
        pD["size"] = "-1"
        pD["accession"] = ",".join(idList)
        ureq = UrlRequestUtil()
        logger.info("doRequestSecondary:  %r  %r  %r  %r", baseUrl, endPoint, pD, hL)
        return ureq.get(baseUrl, endPoint, pD, headers=hL)

    def __doMappingRequest(self, idList):
        """ """
        baseUrl = self.__urlPrimary
        ureq = UrlRequestUtil()
        endPoint = "idmapping/run"
        pD = {"from": "UniProtKB_AC-ID", "to": "UniProtKB", "ids": ",".join(idList)}
        rspJson, retCode = ureq.post(baseUrl, endPoint, pD, headers=[], returnContentType="JSON")
        if retCode != 200:
            logger.error("Primary request failed with retCode %r", retCode)
            return None, retCode
        jobId = rspJson["jobId"]
        logger.debug("jobId %r", jobId)
        ok = self.__checkIdMappingResultsReady(jobId, timeout=600)
        if not ok:
            logger.error("Job failed to run or never finished.")
            return None, retCode
        endPointResults = os.path.join("idmapping/uniprotkb/results", jobId)
        hD = {"Accept": "application/xml"}
        ret, respCode = ureq.getUnWrapped(baseUrl, endPointResults, paramD=None, headers=hD, overwriteUserAgent=False)

        return ret, respCode

    def __processIdList(self, idList):
        """
        Filter the input id list for variants and create the list of unique
        searchable id codes.   Create a dictionary of variant identifiers and
        their corresponding searchable id codes.

        """
        variantD = {}
        tList = []
        #
        for tId in idList:
            # check for variant id
            idx = tId.find("-")
            if idx == -1:
                sId = tId
            else:
                sId = tId[0:idx]
                variantD[tId] = sId
            #
            tList.append(sId)
        #
        #
        return tList, variantD

    def __makeSubLists(self, num, iterable):
        args = [iter(iterable)] * num
        return ([e for e in t if e is not None] for t in zip_longest(*args))

    def __makeSubListsWithPadding(self, num, iterable, padvalue=None):
        "__sublist(3, 'abcdefg', 'x') --> ('a','b','c'), ('d','e','f'), ('g','x','x')"
        return zip_longest(*[iter(iterable)] * num, fillvalue=padvalue)

    def fetchSequenceList(self, unpIdList, retryAltApi=True, usePrimary=True):
        return self.__doSequenceRequest(unpIdList, retryAltApi=retryAltApi, usePrimary=usePrimary)

    def __doSequenceRequest(self, unpIdList, retryAltApi=True, usePrimary=True):
        ok = False
        sD = {}
        if usePrimary:
            ok, sD = self.__doSequenceRequestPrimary(unpIdList)
        #
        if retryAltApi and not ok:
            tIdList = set(unpIdList) - set(sD.keys())
            logger.info("Retrying using secondary service site for (%d) id codes", len(tIdList))
            ok, tD = self.__doSequenceRequestSecondary(tIdList)
            if tD:
                sD.update(tD)
        #
        #
        return ok, sD

    def __doSequenceRequestPrimary(self, unpIdList):
        """ """
        sD = {}
        baseUrl = self.__urlPrimary
        hD = {"Accept": "text/x-fasta"}
        pD = {"format": "fasta"}
        ureq = UrlRequestUtil()
        ok = True
        for unpId in unpIdList:
            endPoint = "uniprotkb/" + unpId
            ret, retCode = ureq.getUnWrapped(baseUrl, endPoint, pD, headers=hD)
            logger.debug("unpId %r url %s endpoint %r ret %r retCode %r", unpId, baseUrl, endPoint, ret, retCode)
            if retCode in [200] and ret and len(ret) > 0:
                rOk, seqId, rD = self.__parseFastaResponse(ret)
                if rOk:
                    sD[seqId] = rD
                else:
                    logger.error("Parsing error in sequence data for %r", unpId)
            else:
                ok = False
        return ok, sD

    def __doSequenceRequestSecondary(self, unpIdList):
        """ """
        sD = {}
        baseUrl = self.__urlSecondary
        hD = {"Accept": "text/x-fasta"}
        pD = {}
        ok = True
        for unpId in unpIdList:
            endPoint = "proteins/api/proteins/" + unpId
            ureq = UrlRequestUtil()
            ret, retCode = ureq.getUnWrapped(baseUrl, endPoint, pD, headers=hD)
            if retCode in [200] and ret and len(ret) > 0:
                rOk, seqId, rD = self.__parseFastaResponse(ret)
                if rOk:
                    sD[seqId] = rD
                else:
                    logger.error("Parsing error in sequence data for %r", unpId)
            else:
                ok = False
        return ok, sD

    def __parseFastaResponse(self, rspS):
        """Parse the text response for FASTA data -

        Args:
            rspS (str): string response containing FASTA data

        Retuns:
            (seqId, dict):  {seqId: {'sequence': 'AAAAAA', 'ky1": va1, ... }}
        """
        faU = FastaUtil()
        cD = {}
        try:
            rL = rspS.splitlines()
            seqId, cD = faU.parseComment(rL[0], "uniprot_regex")
            ok, cD["sequence"] = faU.cleanSequence("".join(rL[1:]))
            return ok, seqId, cD
        except Exception as e:
            logger.exception("Failing with %s", str(e))
        return None, cD

    def doLookup(self, itemList, itemKey="Gen_Name"):
        """Do accession code lookup by mapping from itemKey to accession code for itemList.
        """
        rL = []
        try:
            baseUrl = self.__urlPrimary
            endPoint = "idmapping/run"
            hL = []
            pD = {"from": itemKey, "to": "UniProtKB", "ids": ",".join(itemList)}
            ureq = UrlRequestUtil()
            rspJson, retCode = ureq.post(baseUrl, endPoint, pD, headers=hL, returnContentType="JSON")
            jobId = rspJson["jobId"]
            ok = self.__checkIdMappingResultsReady(jobId, timeout=600)
            if not ok:
                logger.error("Job failed to run or never finished.")
            else:
                batchUrl = os.path.join(self.__urlPrimary, "idmapping/results", jobId)
                rL = []
                for batch in self.__getBatch(batchUrl):
                    resD = batch.json()["results"]
                    idList = [rD["to"] for rD in resD]
                    rL += idList
            return rL, retCode
        except Exception as e:
            logger.exception("Failing with %s", str(e))
        return rL, None

    def __checkIdMappingResultsReady(self, jobId, checkInterval=1, timeout=600):
        """Check status of UniProt REST request job.

        Args:
            jobId (str): Job ID
            checkInterval (int, optional): number seconds to wait between request attempts. Defaults to 10.
            timeout (int, optional): max seconds to wait for job to complete. Defaults to 600.
        """
        try:
            baseUrl = self.__urlPrimary
            endPoint = os.path.join("idmapping/status", jobId)
            ureq = UrlRequestUtil()
            waitTime = 0
            while waitTime < timeout:
                rspJson, retCode = ureq.get(baseUrl, endPoint, {}, headers=[], returnContentType="JSON")
                if retCode != 200:
                    logger.error("Request failed with return code %r and response %r", retCode, rspJson)
                if "jobStatus" in rspJson:
                    if rspJson["jobStatus"] == "RUNNING":
                        logger.info("Retrying in %r s", checkInterval)
                        time.sleep(checkInterval)
                    else:
                        raise Exception(rspJson["jobStatus"])
                else:
                    response = requests.get(os.path.join(baseUrl, endPoint), timeout=600)
                    if "results" in rspJson:
                        total = response.headers["X-Total-Results"]
                        logger.info("Total number of resulting IDs mapped: %r", total)
                    if "failedIds" in rspJson:
                        logger.info("Total number of failedIds: %r", len(rspJson["failedIds"]))
                    return bool(rspJson["results"] or rspJson["failedIds"])
        except Exception as e:
            logger.exception("Failing with %s", str(e))
        return False

    def __getNextLink(self, headers):
        """Get link for the next pagination of results.

        Code from: https://www.uniprot.org/help/api_queries#12-large-number-of-results-use-pagination

        Args:
            headers (dict): response headers
        """
        reNextLink = re.compile(r'<(.+)>; rel="next"')
        if "Link" in headers:
            match = reNextLink.match(headers["Link"])
            if match:
                return match.group(1)

    def __getBatch(self, batchUrl):
        while batchUrl:
            logger.debug("batchUrl %r:", batchUrl)
            response = requests.get(batchUrl, timeout=600)
            response.raise_for_status()
            respCode = response.status_code
            if respCode != 200:
                logger.error("Batch request returned non-200 response code: %r", respCode)
            yield response
            batchUrl = self.__getNextLink(response.headers)

    def doGeneLookup(self, geneName, taxId, reviewed=False):
        """ """
        rL = []
        try:
            baseUrl = self.__urlPrimary
            endPoint = "uniprotkb/search"
            hL = []
            if reviewed:
                pD = {"query": 'gene:"%s" AND taxonomy_id:%s AND reviewed:yes' % (geneName, taxId), "format": "list"}
            else:
                pD = {"query": 'gene:"%s" AND taxonomy_id:%s' % (geneName, taxId), "format": "list"}
            ureq = UrlRequestUtil()
            rspTxt, retCode = ureq.get(baseUrl, endPoint, pD, headers=hL)
            tValL = rspTxt.split("\n") if rspTxt else []
            idList = [tVal for tVal in tValL if tVal]
            return idList, retCode
        except Exception as e:
            logger.exception("Failing with %s", str(e))
        return rL, None
