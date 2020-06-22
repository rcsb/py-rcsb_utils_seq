##
# File:    SeqAlign.py
# Author:  J. Westbrook
# Date:    20-Sep-2019
#
# Updates:

##
"""
Utilities for processing sequence alignments.

"""
__docformat__ = "restructuredtext en"
__author__ = "John Westbrook"
__email__ = "jwest@rcsb.rutgers.edu"
__license__ = "Apache 2.0"

import logging
import sys

if sys.version_info[0] < 3:
    from backports.range import range  # pylint: disable=no-name-in-module,import-error,redefined-builtin

logger = logging.getLogger(__name__)


def doRangesOverlap(r1, r2):
    if r1.start == r1.stop or r2.start == r2.stop:
        return False
    return (r1.start < r2.stop and r1.stop > r2.start) or (r1.stop > r2.start and r2.stop > r1.start)


def getRangeOverlap(r1, r2):
    if not doRangesOverlap(r1, r2):
        return set()
    return set(range(max(r1.start, r2.start), min(r1.stop, r2.stop) + 1))


def splitSeqAlignObjList(seqAlignObjL):
    """Separate the input list range objects into sublists of non-overlapping range segments

    Args:
        saObjL (list): list of SequenceAlignment objects

    Returns:
        (dict): dictionary of sublists (w/ keys 1,2,3) of non-overlapping SequenceAlignment objects
    """
    grpD = {}
    numG = 0
    try:
        seqAlignObjL.sort(key=lambda saObj: saObj.getEntityRange().stop - saObj.getEntityRange().start, reverse=True)
        for saObj in seqAlignObjL:
            inGroup = False
            igrp = 0
            for grp, tsaObjL in grpD.items():
                inGroup = any([doRangesOverlap(saObj.getEntityRange(), tsaObj.getEntityRange()) for tsaObj in tsaObjL])
                if inGroup:
                    igrp = grp
                    break
            numG = numG if inGroup else numG + 1
            igrp = igrp if inGroup else numG
            grpD.setdefault(igrp, []).append(saObj)
    except Exception as e:
        logger.exception("Failing with %s", str(e))
    return grpD


class SeqAlign(object):
    """

    """

    def __init__(self, alignType, **kwargs):
        if alignType == "PDB":
            self.__entitySeqIdBeg = kwargs.get("entitySeqIdBeg", None)
            self.__entitySeqIdEnd = kwargs.get("entitySeqIdEnd", None)
            self.__entityAlignLength = kwargs.get("entityAlignLength", None)
            self.__dbSeqIdBeg = kwargs.get("dbSeqIdBeg", None)
            self.__dbSeqIdEnd = kwargs.get("dbSeqIdEnd", None)
            self.__dbName = kwargs.get("dbName", None)
            self.__dbAccession = kwargs.get("dbAccession", None)
            self.__dbIsoform = kwargs.get("dbIsoform", None)
        elif alignType == "SIFTS":
            # {"UP": unpId, "BG": entitySeqBeg, "UBG": unpSeqBeg, "LEN": entityLength}
            self.__entitySeqIdBeg = kwargs.get("BG", None)
            self.__entityAlignLength = kwargs.get("LEN", None)
            self.__entitySeqIdEnd = self.__entitySeqIdBeg + self.__entityAlignLength - 1
            self.__dbSeqIdBeg = kwargs.get("UBG", None)
            self.__dbSeqIdEnd = kwargs.get("UND", None)
            self.__dbName = "UNP"
            self.__dbAccession = kwargs.get("UP", None)
            self.__dbIsoform = None
        #

    def getEntityRange(self):
        try:
            return range(int(self.__entitySeqIdBeg), int(self.__entitySeqIdEnd) + 1)
        except Exception:
            pass
        return None

    def getEntitySeqIdBeg(self):
        try:
            return int(self.__entitySeqIdBeg)
        except Exception:
            pass
        return None

    def getEntitySeqIdEnd(self):
        try:
            return int(self.__entitySeqIdEnd)
        except Exception:
            pass
        return None

    def getDbSeqIdBeg(self):
        try:
            return int(self.__dbSeqIdBeg)
        except Exception:
            pass
        return None

    def getDbSeqIdEnd(self):
        try:
            return int(self.__dbSeqIdEnd)
        except Exception:
            pass
        return None

    def getEntityAlignLength(self):
        length = None
        try:
            length = self.__entityAlignLength if self.__entityAlignLength else self.__entitySeqIdEnd - self.__entitySeqIdBeg + 1
        except Exception:
            pass
        return length

    def getDbName(self):
        return self.__dbName

    def getDbAccession(self):
        return self.__dbAccession

    def getDbIsoform(self):
        return self.__dbIsoform

    def __str__(self):
        return "DB: %r ACC: %r ISOFORM %r ENITY BEG: %r DB BEG: %r LEN: %r" % (
            self.__dbName,
            self.__dbAccession,
            self.__dbIsoform,
            self.__entitySeqIdBeg,
            self.__dbSeqIdBeg,
            self.__entityAlignLength,
        )

    def __repr__(self):
        return self.__str__()
