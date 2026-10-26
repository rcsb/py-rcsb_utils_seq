##
#  File:           GlycanProvider.py
#  Date:           28-May-2021 jdw
#
#  Updated:
#  13-Jul-2026 dwp Add getGlycanIdentifier() to look up a glycan identifier by entry ID (short or extended) and entity ID
#  20-Oct-2026 dwp Document that identifier map keys follow the entity ID format of the source database
#                  (replacing outdated extended ID TODOs)
#  26-Oct-2026 dwp Use lazy % formatting in a debug log message (pylint W1203)
##
"""
Accessors for glycan mapped annotations.

"""

import logging
import os.path
import time

from rcsb.utils.io.MarshalUtil import MarshalUtil
from rcsb.utils.io.StashableBase import StashableBase

logger = logging.getLogger(__name__)


class GlycanProvider(StashableBase):
    """Accessors (only) for entity glycan mapped identifiers.

    dirPath -> CACHE/glycan/
                             branched_entity_glycan_identifier_map.json
                             accession-wurcs-mapping.json
                    stash/glycan.tar.gz

        Note: This provider supports the delivery of pre-generated data from
                   rcsb.exdb.branch/GlycanProvider

        The identifier map is keyed by branched entity ID, taken as-is from pdbx_core_branched_entity.rcsb_id in the
        database the data was built from. Keys therefore follow that database's ID format: short IDs in upper case
        (e.g., '1ABC_2') or extended IDs in lower case (e.g., 'pdb_00001abc_2'). No change is needed here or in the
        data prep for the transition to extended IDs.
    """

    def __init__(self, **kwargs):
        #
        self.__version = "0.50"
        cachePath = kwargs.get("cachePath", ".")
        useCache = kwargs.get("useCache", True)
        self.__dirName = "glycan"
        self.__dirPath = os.path.join(cachePath, self.__dirName)
        super(GlycanProvider, self).__init__(cachePath, [self.__dirName])
        #
        self.__mU = MarshalUtil(workPath=self.__dirPath)
        self.__glyD = self.__reload(fmt="json", useCache=useCache)
        #

    def testCache(self, minCount=1):
        if minCount == 0:
            return True
        if self.__glyD and minCount and ("identifiers" in self.__glyD) and len(self.__glyD["identifiers"]) >= minCount:
            logger.info("Glycan identifiers (%d)", len(self.__glyD["identifiers"]))
            return True
        return False

    def getIdentifiers(self):
        """Return a dictionary of related identifiers organized by branched entity id.

        Returns:
            (dict): {entityId: {'idType1': ids, 'idType1': ids}, ... }
        """
        try:
            return self.__glyD["identifiers"] if self.__glyD["identifiers"] else {}
        except Exception as e:
            logger.error("Failing with %r", str(e))
        return {}

    def getGlycanIdentifier(self, entryId, entityIdSuffix, idTypeFilter):
        """Get the Glycan identifier for the given ID type (e.g., "glyTouCanId").

        Args:
            entryId (str): entry ID, either short (4 character) or extended (12 character). Lookups are case-sensitive,
                           so pass it in the case stored in the identifier map: upper case for short IDs (e.g., '1ABC'),
                           lower case for extended IDs (e.g., 'pdb_00001abc')
            entityIdSuffix (str): entity ID suffix (e.g., "1", "2", ...)
            idTypeFilter (str): type of ID to return (e.g., "glyTouCanId")

        Returns:
            str: Gly ID of the provided ID type, or None if not found
        """
        gId = None
        glyIdD = self.getIdentifiers()
        try:
            entityId = f"{entryId}_{entityIdSuffix}"
            gIdD = glyIdD.get(entityId, None)
            if gIdD and idTypeFilter in gIdD:
                gId = gIdD.get(idTypeFilter, None)
            if gIdD is not None:
                logger.debug("Got gId %s from entity ID %s", gId, entityId)
        except Exception:
            pass
        return gId

    def __getMappingFilePath(self, fmt="json"):
        baseFileName = "branched_entity_glycan_identifier_map"
        fExt = ".json" if fmt == "json" else ".pic"
        fp = os.path.join(self.__dirPath, baseFileName + fExt)
        return fp

    def reload(self):
        """Reload from the current cache file."""
        ok = False
        try:
            self.__glyD = self.__reload(fmt="json", useCache=True)
            ok = self.__glyD is not None
        except Exception as e:
            logger.exception("Failing with %s", str(e))
        return ok

    def __reload(self, fmt="json", useCache=True):
        mappingFilePath = self.__getMappingFilePath(fmt=fmt)
        tS = time.strftime("%Y %m %d %H:%M:%S", time.localtime())
        pcD = {"version": self.__version, "created": tS, "identifiers": {}}

        if useCache and self.__mU.exists(mappingFilePath):
            logger.info("reading cached path %r", mappingFilePath)
            pcD = self.__mU.doImport(mappingFilePath, fmt=fmt)
        return pcD
