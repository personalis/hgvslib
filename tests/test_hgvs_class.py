#!/bin/py
import unittest


from hgvslib.class_functions import check_hgvs_status, compare_hgvs
from hgvslib.effect import Effect
from hgvslib import constants as c


#----------------------------------------------------
# Test discimination between coding and protein hgvs
#-------------------------------------------------------
class test_compare_hgvs(unittest.TestCase):

	def setUp(self):
		self.longMessage=True

	def test_compare_hgvs_1(self):
		self.hgvs1='c.652delG'
		self.hgvs2='c.652del'
		self.check = c.EQUIVALENT
		check = compare_hgvs(self.hgvs1, self.hgvs2)
		print('Check %s..' % (self._testMethodName))
		self.assertEqual(self.check, check )


	def test_compare_hgvs_2(self):
		self.hgvs1='p.Arg652*'
		self.hgvs2='p.Arg652Ter'
		self.check = c.EQUIVALENT

		check = compare_hgvs(self.hgvs1, self.hgvs2)
		print('Check %s..' % (self._testMethodName))
		self.assertEqual(self.check, check )

	def test_compare_hgvs_3(self):
		self.hgvs1='p.Ala859Ala='
		self.hgvs2='p.Ala859Ala'
		self.check = c.EQUIVALENT

		check = compare_hgvs(self.hgvs1, self.hgvs2)
		print('Check %s..' % (self._testMethodName))
		self.assertEqual(self.check, check )

	def test_compare_hgvs_4(self):
		self.hgvs1='p.Ala859='
		self.hgvs2='p.='
		self.check = c.EQUIVALENT

		check = compare_hgvs(self.hgvs1, self.hgvs2)
		print('Check %s..' % (self._testMethodName))
		self.assertEqual(self.check, check )

	def test_compare_hgvs_5(self):
		self.hgvs1='p.Ala859Ala'
		self.hgvs2='p.(=)'
		self.check = c.EQUIVALENT

		check = compare_hgvs(self.hgvs1, self.hgvs2)
		print('Check %s..' % (self._testMethodName))
		self.assertEqual(self.check, check )

	def test_compare_hgvs_6(self):
		self.hgvs1='p.Gln921Ter'
		self.hgvs2='p.Gln921*'
		self.check = c.EQUIVALENT

		check = compare_hgvs(self.hgvs1, self.hgvs2)
		print('Check %s..' % (self._testMethodName))
		self.assertEqual(self.check, check )

	def test_compare_hgvs_7(self):
		self.hgvs1='p.Ter257Ter'
		self.hgvs2='p.Ter257='
		self.check = c.EQUIVALENT

		check = compare_hgvs(self.hgvs1, self.hgvs2)
		print('Check %s..' % (self._testMethodName))
		self.assertEqual(self.check, check )

	def test_compare_hgvs_8(self):
		self.hgvs1='p.Ter257Ter'
		self.hgvs2='p.*257*'
		self.check = c.EQUIVALENT

		check = compare_hgvs(self.hgvs1, self.hgvs2)
		print('Check %s..' % (self._testMethodName))
		self.assertEqual(self.check, check )


#----------------------------------------------------------------------------------------------------------------------
# check compare protein and coding hgvs function
#----------------------------------------------------------------------------------------------------------------------

class test_compare_protein_coding_hgvs(unittest.TestCase):


	def setUp(self):
		self.longMessage=True

	#--------------------------------------------
	# check protein frameshift
	#--------------------------------------------

	def test_compare_function_on_frameshift_1(self):
		self.query = 'p.Glu78fs'
		self.ref= 'p.Glu78Glyfs*7'
		self.check = c.EQUIVALENT

		message = '%s, %s, %s, %s.' % (self._testMethodName, self.query, self.ref, self.check)
		self.assertEqual(compare_hgvs(self.ref, self.query), self.check)

	def test_compare_function_on_frameshift_2(self):
		self.query ='p.Glu78Gly'
		self.ref='p.Glu78Glyfs*7'
		self.check='no'

		message = '%s, %s, %s, %s.' % (self._testMethodName, self.query, self.ref, self.check)
		self.assertEqual(compare_hgvs(self.query, self.ref ), self.check)

	def test_compare_function_on_frameshift_3(self):
		self.query ='p.Gln188fs'
		self.ref='p.Gln188Glufs'
		self.check = c.EQUIVALENT
		message = '%s, %s, %s, %s.' % (self._testMethodName, self.query, self.ref, self.check)
		self.assertEqual(compare_hgvs(self.query, self.ref ), self.check)


	def test_compare_function_on_frameshift_4(self):
		self.query ='p.Arg1921fs'
		self.ref='-'
		self.check='no'

		message = '%s, %s, %s, %s.' % (self._testMethodName, self.query, self.ref, self.check)
		self.assertEqual(compare_hgvs(self.query, self.ref ), self.check)

	def test_compare_function_on_frameshift_5(self):
		self.query ='p.Thr655Asnfs'
		self.ref='p.Thr655Asnfs*49'
		self.check = c.EQUIVALENT

		message = '%s, %s, %s, %s.' % (self._testMethodName, self.query, self.ref, self.check)
		self.assertEqual(compare_hgvs(self.query, self.ref ), self.check)

	def test_compare_function_on_frameshift_6(self):
		self.query ='p.Arg1921Metfs*'
		self.ref='p.Arg1921fs'
		self.check = c.EQUIVALENT

		message = '%s, %s, %s, %s.' % (self._testMethodName, self.query, self.ref, self.check)
		self.assertEqual(compare_hgvs(self.query, self.ref ), self.check)

	def test_compare_function_on_frameshift_7(self):
		self.query ='p.Arg1921Metfs*'
		self.ref='p.Arg1921MetfsTer9'
		self.check = c.EQUIVALENT

		message = '%s, %s, %s, %s.' % (self._testMethodName, self.query, self.ref, self.check)
		self.assertEqual(compare_hgvs(self.query, self.ref ), self.check)

	def test_compare_function_on_frameshift_8(self):
		self.query ='p.Asn1571fs'
		self.ref='p.Asn1571fs*40'
		self.check = c.EQUIVALENT

		message = '%s, %s, %s, %s.' % (self._testMethodName, self.query, self.ref, self.check)
		self.assertEqual(compare_hgvs(self.query, self.ref ), self.check)



	def test_compare_function_on_frameshift_9(self):
		self.query ='p.Gln64fs*>182'
		self.ref='p.Gln64LysfsTer218'
		self.check = c.EQUIVALENT

		message = '%s, %s, %s, %s.' % (self._testMethodName, self.query, self.ref, self.check)
		self.assertEqual(compare_hgvs(self.query, self.ref ), self.check)


	def test_compare_function_on_nonsense_7(self):
		self.query ='p.Arg168Ter'
		self.ref='p.Arg168*'
		self.check = c.EQUIVALENT

		message = '%s, %s, %s, %s.' % (self._testMethodName, self.query, self.ref, self.check)
		self.assertEqual(compare_hgvs(self.query, self.ref ), self.check)

	def test_compare_function_on_ext_1(self):
		self.query ='p.Ter330Trpext*?'
		self.ref='p.Ter330TrpextTer19'
		self.check = c.EQUIVALENT
		message = '%s, %s, %s, %s.' % (self._testMethodName, self.query, self.ref, self.check)
		self.assertEqual(compare_hgvs(self.query, self.ref ), self.check)

	def test_compare_function_on_ext_2(self):
		self.query ='p.Ter330Trp'
		self.ref='p.Ter330Trpext*?'
		self.check = c.EQUIVALENT
		message = '%s, %s, %s, %s.' % (self._testMethodName, self.query, self.ref, self.check)
		self.assertEqual(compare_hgvs(self.query, self.ref ), self.check)

	def test_compare_function_on_ext_3(self):
		self.query ='p.Ter330TrpextTer19'
		self.ref='p.Ter330Trp'
		self.check = c.EQUIVALENT
		message = '%s, %s, %s, %s.' % (self._testMethodName, self.query, self.ref, self.check)
		self.assertEqual(compare_hgvs(self.query, self.ref ), self.check)

	def test_compare_function_on_synonymous(self):
		self.query ='p.Ter194='
		self.ref='p.Ter194Ter'
		self.check = c.EQUIVALENT
		message = '%s, %s, %s, %s.' % (self._testMethodName, self.query, self.ref, self.check)
		self.assertEqual(compare_hgvs(self.query, self.ref ), self.check)

	def test_compare_function_on_synonymous_2(self):
		self.query ='p.Pro401Pro='
		self.ref='p.Pro401Pro'
		self.check = c.EQUIVALENT
		message = '%s, %s, %s, %s.' % (self._testMethodName, self.query, self.ref, self.check)
		self.assertEqual(compare_hgvs(self.query, self.ref ), self.check)

	def test_compare_function_on_mnv(self):
		self.query ='p.Phe334_Val335delinsLeuMet'
		self.ref='p.PheVal334LeuMet'
		self.check = c.EQUIVALENT
		message = '%s, %s, %s, %s.' % (self._testMethodName, self.query, self.ref, self.check)
		self.assertEqual(compare_hgvs(self.query, self.ref ), self.check)


	#--------------------------------------------------------------------------------
	# check coding deletions
	#--------------------------------------------------------------------------------

	def test_compare_function_on_del_1(self):
		self.query ='c.919_922del'
		self.ref='c.919_922delAG'
		self.check = c.EQUIVALENT

		message = '%s, %s, %s, %s.' % (self._testMethodName, self.query, self.ref, self.check)
		self.assertEqual(compare_hgvs(self.ref, self.query), self.check)


	def test_compare_function_on_del_2(self):
		self.query ='c.302+1_302+4delGTGA'
		self.ref='c.302+1_302+4del'
		self.check = c.EQUIVALENT

		message = '%s, %s, %s, %s.' % (self._testMethodName, self.query, self.ref, self.check)
		self.assertEqual(compare_hgvs(self.ref, self.query), self.check)


	# test preferred and non preferred deletion hgvs equivalencies
	def test_compare_function_on_del_3(self):
		self.query='c.843_846+4delCGAGGTGA'
		self.ref ='c.843_846+4del8'
		self.check = c.EQUIVALENT

		message = '%s, %s, %s, %s.' % (self._testMethodName, self.query, self.ref, self.check)
		self.assertEqual(compare_hgvs(self.ref, self.query), self.check)

	def test_compare_function_on_del_4(self):
		self.query='c.843_846+4delCGAGGTGA'
		self.ref ='c.843_846+4del8'
		self.check = c.EQUIVALENT

		message = '%s, %s, %s, %s.' % (self._testMethodName, self.query, self.ref, self.check)
		self.assertEqual(compare_hgvs(self.ref, self.query), self.check)


	def test_compare_function_on_del_5(self):
		self.query='c.209delCinsCA'
		self.ref ='c.210delA'
		self.check='no'

		message = '%s, %s, %s, %s.' % (self._testMethodName, self.query, self.ref, self.check)
		self.assertEqual(compare_hgvs(self.ref, self.query), self.check)

	def test_compare_function_on_del_6(self):
		self.query='c.3324_3347delTGCAGCTGCAGCTGCAGCCGCAGCTGC'
		self.ref ='c.3324_3347delAGCTGCTGCAGCTGCAGCTGCAGC'
		self.check='no'

		message = '%s, %s, %s, %s.' % (self._testMethodName, self.query, self.ref, self.check)
		self.assertEqual(compare_hgvs(self.ref, self.query), self.check)


	def test_compare_function_on_ins_7(self):
		self.query='NM_033632.3:c.45_46insCCT'
		self.ref ='c.45_46insCCT'
		self.check = c.EXACT

		message = '%s, %s, %s, %s.' % (self._testMethodName, self.query, self.ref, self.check)
		self.assertEqual(compare_hgvs(self.ref, self.query), self.check)


#--------------------------------------------
# check effect impact
#--------------------------------------------

class test_check_effect_impact(unittest.TestCase):


	def setUp(self):
		self.longMessage=True

	def test_effect_missense_variant_1(self):
		self.query = 'missense_variant'
		self.ref='intron_variant'
		self.check='no'

		message = '%s, %s, %s, %s.' % (self._testMethodName, self.query, self.ref, self.check)
		self.assertEqual(Effect.check_effect(self.query, self.ref), self.check)

	def test_effect_missense_2(self):
		query = 'missense_variant'
		ref = 'non_synonymous_codon'
		check = c.EQUIVALENT
		message = '%s, %s, %s, %s.' % (self._testMethodName, query, ref, check)
		self.assertEqual(Effect.check_effect(query, ref),check)

	def test_effect_insertion_1(self):
		query = 'inframe_insertion'
		ref = 'disruptive_inframe_insertion'
		check = c.EQUIVALENT
		message = '%s, %s, %s, %s.' % (self._testMethodName, query, ref, check)
		self.assertEqual(Effect.check_effect(query, ref),check)

	def test_effect_insertion_2(self):
		query = 'inframe_insertion'
		ref = 'inframe_insertion'
		check = c.EXACT
		message = '%s, %s, %s, %s.' % (self._testMethodName, query, ref, check)
		self.assertEqual(Effect.check_effect(query, ref),check)
