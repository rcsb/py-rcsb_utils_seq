##
# File:    SiftsUtils.py
# Author:  J. Westbrook
# Date:    11-Dec-2018
#
# Updates:

##
"""
Utilities to read SIFTS PDB -> UniProt mapping data.

"""
__docformat__ = "restructuredtext en"
__author__ = "John Westbrook"
__email__ = "jwest@rcsb.rutgers.edu"
__license__ = "Apache 2.0"

import logging
import os

from rcsb.utils.io.MarshalUtil import MarshalUtil

logger = logging.getLogger(__name__)


class SiftsUtils(object):
    """
    """

    def __init__(self):
        self.__ns = None

    def __readSiftsSummaryFile(self, filePath):
        """ Read input SIFTS summary file and return a list of dictionaries.
        """
        cL = []
        try:
            mU = MarshalUtil()
            cL = mU.doImport(filePath, format="csv")
            logger.debug("Container list %d" % len(cL))
        except Exception as e:
            logger.exception("Failing with %s" % str(e))
        return cL

    def getSummaryMapping(self, siftsSummaryDirPath):
        """ Integrated SIFTS summary instance-level taxonomy and UniProt mapping data.

        Returns:


        # Taxonomy
        [('PDB', '101m'), ('CHAIN', 'A'), ('TAX_ID', '9755'), ('SCIENTIFIC_NAME', 'PHYCD')]
        # Uniprot
        [('PDB', '101m'), ('CHAIN', 'A'), ('SP_PRIMARY', 'P02185'), ('RES_BEG', '1'), ('RES_END', '154'), ('PDB_BEG', '0'), ('PDB_END', '153'), ('SP_BEG', '1'), ('SP_END', '154')]
        """
        fp = os.path.join(siftsSummaryDirPath, 'pdb_chain_taxonomy.csv.gz')
        taxDL = self.__readSiftsSummaryFile(fp)
        logger.info("Length of SIFTS taxonomy summary file %d" % len(taxDL))
        logger.info("%r" % list(taxDL[0].items()))
        tD = {}
        for taxD in taxDL:
            entryId = taxD['PDB']
            chainId = taxD['CHAIN']
            taxId = taxD['TAX_ID']
            #
            if entryId not in tD:
                tD[entryId] = {}
            #
            if chainId not in tD[entryId]:
                tD[entryId][chainId] = {}
            #
            tD[entryId][chainId][taxId] = True
        #
        logger.info("Taxonomy for %d entries" % len(tD))
        #
        #
        fp = os.path.join(siftsSummaryDirPath, 'pdb_chain_uniprot.csv.gz')
        unpDL = self.__readSiftsSummaryFile(fp)
        logger.info("Length of SIFTS UniProt summary file %d" % len(unpDL))
        logger.info("%r" % list(unpDL[0].items()))
        uD = {}
        for unpD in unpDL:
            entryId = unpD['PDB']
            chainId = unpD['CHAIN']
            unpId = unpD['SP_PRIMARY']
            #
            txL = []
            try:
                txL = list(tD[entryId][chainId].keys())
            except Exception:
                logger.debug("No taxonmy for %s %s" % (entryId, chainId))

            #
            pdbSeqBeg = int(unpD['RES_BEG']) if unpD['RES_BEG'].isdigit() else None
            pdbSeqEnd = int(unpD['RES_END']) if unpD['RES_END'].isdigit() else None
            unpSeqBeg = int(unpD['SP_BEG']) if unpD['SP_BEG'].isdigit() else None
            unpSeqEnd = int(unpD['SP_END']) if unpD['SP_END'].isdigit() else None
            #
            if entryId not in uD:
                uD[entryId] = {}
            #
            if chainId not in uD[entryId]:
                uD[entryId][chainId] = []
            #
            uD[entryId][chainId].append({'UP': unpId, 'TX': txL, 'NT': len(txL), 'BG': pdbSeqBeg, 'ND': pdbSeqEnd, 'UBG': unpSeqBeg, 'UND': unpSeqEnd})
            # uD[entryId][chainId].append({'UP': unpId, 'TX': txL, 'NT': len(txL), 'BEG': pdbSeqBeg, 'END': pdbSeqEnd})
            #
        logger.info("UniProt mapping for %d entries" % len(uD))
        # -----
        return uD
