#!/bin/py
import unittest
import re
import sys, os
import logging

from hgvslib.transcript import Transcript


class test_transcript_object(unittest.TestCase):

	def setUp(self):
		self.longMessage=True
		
		self.transcript ="NM_000546.5"
		self.version = '5'
		self.refseq = 'NM_000546'

	def test_transcript_version(self):
		trans = Transcript(self.transcript)

		print("Check %s.." % (self._testMethodName))
		self.assertEqual((self.refseq, self.version), (trans.accession, trans.version))