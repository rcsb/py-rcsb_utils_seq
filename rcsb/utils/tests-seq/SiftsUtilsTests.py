##
# File:    SiftsUtilsTests.py
# Author:  J. Westbrook
# Date:    11-Dec-2018
#
# Update:
#
#
##
"""
Various utilities for processing SIFTS UniProt and taxonomy correspondence data.

"""

__docformat__ = "restructuredtext en"
__author__ = "John Westbrook"
__email__ = "jwest@rcsb.rutgers.edu"
__license__ = "Apache 2.0"

import json
import logging
import os
import unittest

from rcsb.utils.seq.SiftsUtils import SiftsUtils

HERE = os.path.abspath(os.path.dirname(__file__))
TOPDIR = os.path.dirname(os.path.dirname(HERE))

logging.basicConfig(level=logging.INFO)
logger = logging.getLogger()


class SiftsUtilsTests(unittest.TestCase):

    def setUp(self):
        self.__dirPath = os.path.join(os.path.dirname(TOPDIR), 'rcsb', 'mock-data')
        self.__siftsSummaryPath = os.path.join(self.__dirPath, 'sifts-summary')
        #
        self.__siftsMappingFile = os.path.join(HERE, 'test-output', 'uniprot_taxonomy_mapping.json')

    def tearDown(self):
        pass

    def testReadSiftsSummary(self):
        su = SiftsUtils()
        rD = su.getSummaryMapping(self.__siftsSummaryPath)
        #
        logger.info("Model match length %d" % len(rD))

        self.__serializeJson(self.__siftsMappingFile, rD)
        #
        #

    def __serializeJson(self, filePath, oD):
        with open(filePath, "w") as outfile:
            json.dump(oD, outfile, indent=0)


def readSiftsInfo():
    suiteSelect = unittest.TestSuite()
    suiteSelect.addTest(SiftsUtilsTests("testReadSiftsSummary"))
    return suiteSelect


if __name__ == '__main__':

    if True:
        mySuite = readSiftsInfo()
        unittest.TextTestRunner(verbosity=2).run(mySuite)
