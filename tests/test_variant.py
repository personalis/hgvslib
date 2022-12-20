#!/bin/py
import unittest
import re
import logging

from hgvslib.variant import Variant

#--------------------------------------------
# test Variant, Transcript and HGVS creation
#--------------------------------------------

class test_variant_object(unittest.TestCase):

	def setUp(self):
		self.longMessage=True

		self.chgvs ="NM_000546.5:c.652_654delGTG"
		self.phgvs="NP_000537.3:p.Pro559Ser"
		self.id = 'XX'
		self.transcript = 'NM_000546.5'
		self.c_hgvs = 'c.652_654delGTG'
		self.p_hgvs = 'p.Pro559Ser'

	def test_variant(self):
		var = Variant(self.chgvs, self.phgvs, self.id, self.transcript)

		print("Check %s.." % (self._testMethodName))
		self.assertEqual((self.transcript, self.c_hgvs, self.p_hgvs), (var.transcript, var.chgvs, var.phgvs))


class test_from_VCF_info_fields(unittest.TestCase):

	def setUp(self):
		self.longMessage=True 
		
		self.transcript='ENST00000372548'
		self.c_hgvs='c.997C>T'
		self.p_hgvs='p.Arg333*'
		self.protein=''
		self.alt='T'

	def test_variant(self):
		snpeff_info_field = "T|stop_gained|HIGH|CXorf57|ENSG00000147231|transcript|ENST00000372548|protein_coding|4/14|c.997C>T|p.Arg333*|1106/3861|997/2568|333/855||"
		id = 'X1'
		var = Variant.from_vcf_info_field(id, snpeff_info_field)
		print("Check %s.." % (self._testMethodName))
		self.assertEqual((self.transcript, self.c_hgvs, self.p_hgvs, self.alt), (var.transcript, var.chgvs, var.phgvs, var.alt))



