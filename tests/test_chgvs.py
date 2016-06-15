#!/bin/py
import unittest
import re
import sys, os
import logging

from hgvslib.cHGVS import cHGVS
from hgvslib import constants as c


#-------------------------------------------------------------------------------------------------
#  Check coding HGVS functions
#-------------------------------------------------------------------------------------------------


class test_check_cHGVS_type(unittest.TestCase):

	def setUp(self):
		self.longMessage=True
		self.query ='c.22delC'
		self.check='del'

	def test_class_instance_function(self):
		hgvs1=cHGVS(self.query)
		hgvs2=cHGVS(self.query)

		self.assertEqual(hgvs1.name, hgvs2.name)
		self.assertEqual(hgvs1.transcript, hgvs2.transcript)
		self.assertEqual(hgvs1.bases, hgvs2.bases)
		self.assertEqual(hgvs1.alias, hgvs2.alias)
		#self.assertEqual(hgvs1.intronic, hgvs2.intronic)

		self.assertEqual(hgvs1.intronic, hgvs2.intronic)


	def test_class_instance_function_with_transcipt(self):
		query ='NM(GENE):c.22delC'

		hgvs1=cHGVS(query)
		hgvs2=cHGVS(query)

		self.assertEqual(hgvs1.name, hgvs2.name)
		self.assertEqual(hgvs1.transcript, hgvs2.transcript)
		self.assertEqual(hgvs1.bases, hgvs2.bases)
		self.assertEqual(hgvs1.alias, hgvs2.alias)
		self.assertEqual(hgvs1.intronic, hgvs2.intronic)

	def test_get_cHGVS_type_1(self):
		hgvs=cHGVS(self.query)
		message='%s, %s, %s' % (self._testMethodName, self.query, self.check)
		self.assertEqual(hgvs.type, (self.check), msg=message)

	def test_get_cHGVS_type_2(self):
		self.check='delins' #### self.check='delins'
		self.query='c.22_23delTTinsCC'
		hgvs=cHGVS(self.query)
		message='%s, %s, %s' % (self._testMethodName, self.query, self.check)
		self.assertEqual(hgvs.type, (self.check), msg=message)

	def test_get_cHGVS_type_3(self):
		self.check='delins'
		self.query='c.22_23delinsCC'
		hgvs=cHGVS(self.query)
		message='%s, %s, %s' % (self._testMethodName, self.query, self.check)
		self.assertEqual(hgvs.type, (self.check), msg=message)

	def test_get_cHGVS_type_4(self):
		self.query='c.23insCC'
		self.check='ins'
		hgvs=cHGVS(self.query)
		message='%s, %s, %s' % (self._testMethodName, self.query, self.check)
		self.assertEqual(hgvs.type, (self.check), msg=message)

	def test_get_cHGVS_type_5(self):
		self.query='NM_000602.4:c.-820_-817G(4_5)'
		self.check='del'
		hgvs=cHGVS(self.query)
		message='%s, %s, %s' % (self._testMethodName, self.query, self.check)
		self.assertEqual(hgvs.type, (self.check), msg=message)

	def test_get_cHGVS_type_6(self):
		self.query='NM_000240.3:c.-1241_-1212ACCGGCACCGGCACCAGTACCCGCACCAGT(3_5)'
		self.check='del'
		hgvs=cHGVS(self.query)
		message='%s, %s, %s' % (self._testMethodName, self.query, self.check)
		self.assertEqual(hgvs.type, (self.check), msg=message)



#-------------------------------------------------------------------------------------------------
# check delins and dup functions
#-------------------------------------------------------------------------------------------------
class test_check_cHGVS_del_dup_function(unittest.TestCase):

	def setUp(self):
		self.longMessage=True

	def test_check_c_del_dup_1(self):
		self.check = c.EQUIVALENT
		self.query='c.222_223delTT'
		self.ref='c.222_223del'

		hgvs1=cHGVS(self.query)
		hgvs2=cHGVS(self.ref)
		message='%s, %s, %s, %s.' % (self._testMethodName, self.query, self.ref, self.check)
		self.assertEqual(cHGVS.check_c_del_dup(hgvs1, hgvs2 ), self.check)


	def test_check_c_del_dup_2(self):
		self.check = c.EXACT
		self.query='c.222_223delTT'
		self.ref='c.222_223delTT'

		hgvs1=cHGVS(self.query)
		hgvs2=cHGVS(self.ref)
		message='%s, %s, %s, %s.' % (self._testMethodName, self.query, self.ref, self.check)
		self.assertEqual(cHGVS.check_c_del_dup(hgvs1, hgvs2 ), self.check)


	def test_check_c_del_dup_3(self):
		self.check = c.EQUIVALENT
		self.query='c.222_223dup'
		self.ref='c.222_223dupTT'

		hgvs1=cHGVS(self.query)
		hgvs2=cHGVS(self.ref)
		message='%s, %s, %s, %s.' % (self._testMethodName, self.query, self.ref, self.check)
		self.assertEqual(cHGVS.check_c_del_dup(hgvs1, hgvs2 ), self.check)


	def test_c_del_dup_4(self):
		self.check = c.NO_MATCH
		self.query='c.919_922delins'
		self.ref='c.919_922delinsAG'

		hgvs1=cHGVS(self.query)
		hgvs2=cHGVS(self.ref)
		message='%s, %s, %s, %s.' % (self._testMethodName, self.query, self.ref, self.check)
		self.assertEqual(cHGVS.check_c_del_dup(hgvs1, hgvs2 ), self.check)


#-----------------------------------------------------------------------------------------
#  Test the check_c_hgvs function.
#  Note - ins and dup equivalencies are not tested -
#  -- requires going into the fasta file to check that the ins base is a dup
#-----------------------------------------------------------------------------------------

class test_check_cHGVS_comparison(unittest.TestCase):

	def setUp(self):
		self.longMessage=True

	def test_sub_1(self):
		self.check = c.EXACT
		self.query='c.222T>C'
		self.ref='c.222T>C'
		hgvs1=cHGVS(self.query)
		hgvs2=cHGVS(self.ref)
		message='%s, %s, %s, %s.' % (self._testMethodName, self.query, self.ref, self.check)
		self.assertEqual(cHGVS.check_c_hgvs(hgvs1, hgvs2 ), self.check)


	def test_sub_2(self):
		self.check = c.NO_MATCH
		self.query='c.222T>A'
		self.ref='c.222T>C'
		hgvs1=cHGVS(self.query)
		hgvs2=cHGVS(self.ref)
		message='%s, %s, %s, %s.' % (self._testMethodName, self.query, self.ref, self.check)
		self.assertEqual(cHGVS.check_c_hgvs(hgvs1, hgvs2 ), self.check)

	def test_sub_3(self):
		self.check = c.NO_MATCH
		self.query='c.222T>C'
		self.ref='c.222T='
		hgvs1=cHGVS(self.query)
		hgvs2=cHGVS(self.ref)
		message='%s, %s, %s, %s.' % (self._testMethodName, self.query, self.ref, self.check)
		self.assertEqual(cHGVS.check_c_hgvs(hgvs1, hgvs2 ), self.check)



	#--------------------------------------------------------------------------------
	# check indels
	#--------------------------------------------------------------------------------


	def test_delins_1(self):
		self.check = c.EQUIVALENT
		self.query='c.222_223delTTinsGA'
		self.ref='c.222_223delinsGA'
		hgvs1=cHGVS(self.query)
		hgvs2=cHGVS(self.ref)
		message='%s, %s, %s, %s.' % (self._testMethodName, self.query, self.ref, self.check)
		self.assertEqual(cHGVS.check_c_hgvs(hgvs1, hgvs2 ), self.check)


	def test_check_delins_2(self):
		self.check = c.EXACT
		self.query='c.222_223delinsGA'
		self.ref='c.222_223delinsGA'

		hgvs1=cHGVS(self.query)
		hgvs2=cHGVS(self.ref)
		message='%s, %s, %s, %s.' % (self._testMethodName, self.query, self.ref, self.check)
		self.assertEqual(cHGVS.check_c_hgvs(hgvs1, hgvs2 ), self.check)


	def test_check_delins_3(self):
		self.query='c.919_922delins'
		self.ref ='c.919_922delinsAG'
		self.check = c.NO_MATCH

		message='%s, %s, %s, %s.' % (self._testMethodName, self.query, self.ref, self.check)
		var1=cHGVS(self.query)
		var2=cHGVS(self.ref)
		self.assertEqual(cHGVS.check_c_hgvs(var1, var2), self.check)

	def test_check_delins_4(self):
		"""
		Test position of coding hgvs when the actual sequence output does not match.
		"""
		self.query ='c.919_922delinsAC'
		self.ref='c.919_922delinsAG'
		self.check = c.NO_MATCH

		message='%s, %s, %s, %s.' % (self._testMethodName, self.query, self.ref, self.check)
		var1=cHGVS(self.query)
		var2=cHGVS(self.ref)
		self.assertEqual(cHGVS.check_c_hgvs(var1, var2), self.check)


	def test_delins_5(self):
		"""
		Test equivalency of separated del and ins syntax vs delins syntax.
		Preferred hgvs here is describing ins bases only.
		"""
		self.query ='c.711_729delTGAGAGCCGGCTGGCGGATinsCC'
		self.ref='c.711_729delinsCC'
		self.check = c.EQUIVALENT

		message='%s, %s, %s, %s.' % (self._testMethodName, self.query, self.ref, self.check)
		var1=cHGVS(self.query)
		var2=cHGVS(self.ref)
		self.assertEqual(cHGVS.check_c_hgvs(var1, var2), self.check)


	def test_delins_6(self):
		self.query ='c.711_729delTGAGAGCCGGCTGGCGGATinsCCC'
		self.ref='c.711_729delinsCC'
		self.check = c.NO_MATCH

		message='%s, %s, %s, %s.' % (self._testMethodName, self.query, self.ref, self.check)
		self.assertEqual(cHGVS.check_c_hgvs(self.query, self.ref ), self.check)

	def test_delins_7(self):
		ref='c.1077_1093>GCTTTGA'
		query='c.1077_1093delAAinsGCTTTGA'
		check = c.EQUIVALENT
		message='%s, %s, %s, %s.' % (self._testMethodName, ref, query, check)
		var1=cHGVS(query)
		var2=cHGVS(ref)
		self.assertEqual(cHGVS.check_c_hgvs(var1, var2), check)


	def test_delins_8(self):
		ref='c.1077_1093>GCTTTGA'
		query='c.1077_1093delinsGCTTTGA'
		check = c.EQUIVALENT
		message='%s, %s, %s, %s.' % (self._testMethodName, ref, query, check)
		var1=cHGVS(query)
		var2=cHGVS(ref)
		self.assertEqual(cHGVS.check_c_hgvs(var1, var2), check)


	def test_delins_9(self):
		ref='c.1077_1090>GCTTTGA'
		query='c.1077_1093delinsGCTTTGA'
		check='no'
		message='%s, %s, %s, %s.' % (self._testMethodName, ref, query, check)
		var1=cHGVS(query)
		var2=cHGVS(ref)
		self.assertEqual(cHGVS.check_c_hgvs(var1, var2), check)


	#--------------------------------------------------------------------------------
	# check deletions
	#--------------------------------------------------------------------------------

	def test_del_1(self):
		self.query ='c.919_922del'
		self.ref='c.919_922delAG'
		self.check = c.EQUIVALENT

		message='%s, %s, %s, %s.' % (self._testMethodName, self.query, self.ref, self.check)
		var1=cHGVS(self.query)
		var2=cHGVS(self.ref)
		self.assertEqual(cHGVS.check_c_hgvs(var1, var2), self.check)


	def test_del_2(self):
		self.query ='c.302+1_302+4delGTGA'
		self.ref='c.302+1_302+4del'
		self.check = c.EQUIVALENT

		message='%s, %s, %s, %s.' % (self._testMethodName, self.query, self.ref, self.check)
		var1=cHGVS(self.query)
		var2=cHGVS(self.ref)
		self.assertEqual(cHGVS.check_c_hgvs(var1, var2), self.check)



	def test_del_3(self):
		"""
		For preferred and non preferred deletion hgvs equivalencies
		"""
		self.query='c.843_846+4delCGAGGTGA'
		self.ref ='c.843_846+4del8'
		self.check = c.EQUIVALENT

		message='%s, %s, %s, %s.' % (self._testMethodName, self.query, self.ref, self.check)
		var1=cHGVS(self.query)
		var2=cHGVS(self.ref)
		self.assertEqual(cHGVS.check_c_hgvs(var1, var2), self.check)

	def test_del_4(self):
		self.query='c.843_846+4delCGAGGTGA'
		self.ref ='c.843_846+4del8'
		self.check = c.EQUIVALENT

		message='%s, %s, %s, %s.' % (self._testMethodName, self.query, self.ref, self.check)
		var1=cHGVS(self.query)
		var2=cHGVS(self.ref)
		self.assertEqual(cHGVS.check_c_hgvs(var1, var2), self.check)

	def test_del_4(self):
		self.query='c.209delCinsCA'
		self.ref ='c.210delA'
		self.check = c.NO_MATCH

		message='%s, %s, %s, %s.' % (self._testMethodName, self.query, self.ref, self.check)
		var1=cHGVS(self.query)
		var2=cHGVS(self.ref)
		self.assertEqual(cHGVS.check_c_hgvs(var1, var2), self.check)

	def test_del_5(self):
		self.query='c.209_210delCA'
		self.ref ='c.209_210del3'
		self.check = c.NO_MATCH

		message='%s, %s, %s, %s.' % (self._testMethodName, self.query, self.ref, self.check)
		var1=cHGVS(self.query)
		var2=cHGVS(self.ref)
		self.assertEqual(cHGVS.check_c_hgvs(var1, var2), self.check)

	def test_del_6(self):
		self.query='c.209_210delCA'
		self.ref ='c.209_210delCG'
		self.check = c.NO_MATCH
		message='%s, %s, %s, %s.' % (self._testMethodName, self.query, self.ref, self.check)
		var1=cHGVS(self.query)
		var2=cHGVS(self.ref)
		var1=cHGVS(self.query)
		self.assertEqual(cHGVS.check_c_hgvs(var1, var2), self.check)


	#--------------------------------------------------------------------------------
	# check duplications
	#--------------------------------------------------------------------------------


	def test_dup_1(self):
		self.query='c.424-16_424-14dupTTC'
		self.ref='c.424-19_424-14dupTTCTTC'
		self.check = c.NO_MATCH

		message='%s, %s, %s, %s.' % (self._testMethodName, self.query, self.ref, self.check)
		var1=cHGVS(self.query)
		var2=cHGVS(self.ref)
		self.assertEqual(cHGVS.check_c_hgvs(var1, var2), self.check)


	def test_check_dup_2(self):
		"""
		For when the hgvs positions are not even close
		"""
		self.query='c.424-14_424-13insTTC'
		self.ref ='c.424-16_424-14dupTTC'
		self.check = c.NO_MATCH
		message='%s, %s, %s, %s.' % (self._testMethodName, self.query, self.ref, self.check)
		var1=cHGVS(self.query)
		var2=cHGVS(self.ref)
		self.assertEqual(cHGVS.check_c_hgvs(var1, var2), self.check)


	#
	def test_dup_3(self):
		"""
		check dup equivalencies
		"""
		self.query='c.1011dup'
		self.ref='c.1011dupA'
		self.check = c.EQUIVALENT

		message='%s, %s, %s, %s.' % (self._testMethodName, self.query, self.ref, self.check)
		var1=cHGVS(self.query)
		var2=cHGVS(self.ref)
		self.assertEqual(cHGVS.check_c_hgvs(var1, var2), self.check)

	def test_dup_4(self):
		self.query='c.30_38dup'
		self.ref='c.30_38dupGCTGCTGCT'
		self.check = c.EQUIVALENT

		message='%s, %s, %s, %s.' % (self._testMethodName, self.query, self.ref, self.check)
		var1=cHGVS(self.query)
		var2=cHGVS(self.ref)
		self.assertEqual(cHGVS.check_c_hgvs(var1, var2), self.check)



	#--------------------------------------------
	# check insertions
	#--------------------------------------------


	def test_ins_1(self):
		self.query='c.764_765ins33'
		self.ref='c.764_765insAACCTGACAGTTGCAGTTTTCACCCATGGAAAG'
		self.check = c.EQUIVALENT

		message='%s, %s, %s, %s.' % (self._testMethodName, self.query, self.ref, self.check)
		var1=cHGVS(self.query)
		var2=cHGVS(self.ref)
		self.assertEqual(cHGVS.check_c_hgvs(var1, var2), self.check)

	def test_ins_2(self):
		self.query='c.764_765insAACCTGACAGTTGCAGTTTTCACCCATGGAAAGT'
		self.ref='c.764_765insAACCTGACAGTTGCAGTTTTCACCCATGGAAAGT'
		self.check = c.EXACT

		message='%s, %s, %s, %s.' % (self._testMethodName, self.query, self.ref, self.check)
		var1=cHGVS(self.query)
		var2=cHGVS(self.ref)
		self.assertEqual(cHGVS.check_c_hgvs(var1, var2), self.check)


	def test_ins_3(self):
		self.query='c.764_765ins33'
		self.ref='c.764_765ins33'
		self.check = c.EXACT

		message='%s, %s, %s, %s.' % (self._testMethodName, self.query, self.ref, self.check)
		var1=cHGVS(self.query)
		var2=cHGVS(self.ref)
		self.assertEqual(cHGVS.check_c_hgvs(var1, var2), self.check)


class test_cHGVS(unittest.TestCase):

	def setUp(self):
		self.longMessage=True

	def test_cHGVS_parse_hgvs_string(self):
		hgvs_str='NM_017547.3(FOXRED1):c.694C>T'
		var=cHGVS(hgvs_str)
		transcript='NM_017547.3'
		chgvs='c.694C>T'
		print 'Check %s..' % (self._testMethodName)
		self.assertEqual((chgvs, transcript), (var.name, var.transcript))


	def test_is_coding_1(self):
		hgvs_str='NM_017547.3(FOXRED1):c.694C>T'
		print 'Check %s..' % (self._testMethodName)
		self.assertTrue(cHGVS.is_coding(hgvs_str))

	def test_is_coding_2(self):
		hgvs_str='NM_017547.3:m.694C>T'
		print 'Check %s..' % (self._testMethodName)
		self.assertTrue(cHGVS.is_coding(hgvs_str))

	def test_is_coding_3(self):
		hgvs_str='NM_017547.3:r.694C>T'
		print 'Check %s..' % (self._testMethodName)
		self.assertTrue(cHGVS.is_coding(hgvs_str))

	def test_is_coding_4(self):
		hgvs_str='NM_017547.3:n.694C>T'
		print 'Check %s..' % (self._testMethodName)
		self.assertTrue(cHGVS.is_coding(hgvs_str))

	def test_is_coding_5(self):
		hgvs_str='NM_017547.3:g.694C>T'
		print 'Check %s..' % (self._testMethodName)
		self.assertFalse(cHGVS.is_coding(hgvs_str))

	def test_is_coding_6(self):
		hgvs_str='NM_017547.3:p.694C>T'
		print 'Check %s..' % (self._testMethodName)
		self.assertFalse(cHGVS.is_coding(hgvs_str))