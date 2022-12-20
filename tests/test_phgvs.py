#!/bin/py
import unittest
import re
import sys, os, inspect
import logging

from hgvslib.pHGVS import pHGVS


########################################################################################################################
#
# Check protein functions
#
########################################################################################################################


class test_get_protein_variant_type(unittest.TestCase):

	def setUp(self):
		self.longMessage=True

	def test_get_type_1(self):
		self.query ="p.Glu78fs"
		self.check="frameshift"
		var = pHGVS(self.query)
		message = "%s, %s, %s" % (self._testMethodName, self.query, self.check)
		self.assertEqual(var.type, self.check)

	def test_get_type_2(self):
		self.query ="p.Lys2_Met3insGlnSerLys"
		self.check="ins"
		var = pHGVS(self.query)
		message = "%s, %s, %s" % (self._testMethodName, self.query, self.check)
		self.assertEqual(var.type, self.check)


#----------------------------------------------------------------------------------------------------------------------
# check protein hgvs
#----------------------------------------------------------------------------------------------------------------------
class test_check_p_hgvs(unittest.TestCase):

	def setUp(self):
		self.longMessage=True

	#--------------------------------------------
	# check frameshift
	#--------------------------------------------

	def test_frameshift_1(self):
		self.query = "p.Glu78fs"
		self.ref= "p.Glu78Glyfs*7"
		self.refcheck = "p.Glu78fs*"
		self.new_ref = "p.Glu78Glyfs"
		self.check="yes_m"

		hgvs1 = pHGVS(self.query)
		hgvs2 = pHGVS(self.ref)
		message = "%s, %s, %s, %s." % (self._testMethodName, self.query, self.ref, self.check)
		self.assertEqual(pHGVS.check_p_hgvs(hgvs1, hgvs2), self.check)


	def test_frameshift_2(self):
		self.query ="p.Glu78Gly"
		self.ref="p.Glu78Glyfs*7"
		self.check="no"

		hgvs1 = pHGVS(self.query)
		hgvs2 = pHGVS(self.ref)
		message = "%s, %s, %s, %s." % (self._testMethodName, self.query, self.ref, self.check)
		self.assertEqual(pHGVS.check_p_hgvs(self.query, self.ref ), self.check)

	def test_frameshift_3(self):
		self.query ="p.Gln188fs"
		self.ref="p.Gln188Glufs"
		self.check="yes_m"

		hgvs1 = pHGVS(self.query)
		hgvs2 = pHGVS(self.ref)
		message = "%s, %s, %s, %s." % (self._testMethodName, self.query, self.ref, self.check)
		self.assertEqual(pHGVS.check_p_hgvs(self.query, self.ref ), self.check)


	def test_frameshift_4(self):
		self.query ="p.Arg1921fs"
		self.ref="-"
		self.check="no"

		hgvs1 = pHGVS(self.query)
		hgvs2 = pHGVS(self.ref)
		message = "%s, %s, %s, %s." % (self._testMethodName, self.query, self.ref, self.check)
		self.assertEqual(pHGVS.check_p_hgvs(self.query, self.ref ), self.check)

	def test_frameshift_5(self):
		self.query ="p.Thr655Asnfs"
		self.ref="p.Thr655Asnfs*49"
		self.check="yes_m"

		hgvs1 = pHGVS(self.query)
		hgvs2 = pHGVS(self.ref)
		message = "%s, %s, %s, %s." % (self._testMethodName, self.query, self.ref, self.check)
		self.assertEqual(pHGVS.check_p_hgvs(self.query, self.ref ), self.check)

	def test_frameshift_6(self):
		self.query ="p.Arg1921fs"
		self.ref="p.Arg1921MetfsTer9"
		self.check="yes_m"

		hgvs1 = pHGVS(self.query)
		hgvs2 = pHGVS(self.ref)
		message = "%s, %s, %s, %s." % (self._testMethodName, self.query, self.ref, self.check)
		self.assertEqual(pHGVS.check_p_hgvs(self.query, self.ref ), self.check)


	def test_frameshift_6(self):
		self.query ="p.Arg1921fs"
		self.ref="p.Arg1921Metfs"
		self.check="yes_m"

		hgvs1 = pHGVS(self.query)
		hgvs2 = pHGVS(self.ref)
		message = "%s, %s, %s, %s." % (self._testMethodName, self.query, self.ref, self.check)
		self.assertEqual(pHGVS.check_p_hgvs(self.query, self.ref ), self.check)


	#--------------------------------------------
	# check nonsense variants
	#--------------------------------------------

	def test_nonsense_1(self):
		self.query ="p.Tyr567Ter"
		self.ref="p.Tyr567*"
		self.check="yes_m"

		hgvs1 = pHGVS(self.query)
		hgvs2 = pHGVS(self.ref)
		message = "%s, %s, %s, %s." % (self._testMethodName, self.query, self.ref, self.check)
		self.assertEqual(pHGVS.check_p_hgvs(hgvs1, hgvs2 ), self.check)

	def test_nonsense_2(self):
		self.query ="p.Gln136Profs*13"
		self.ref="p.Gln136Ter13"
		self.check="no"

		hgvs1 = pHGVS(self.query)
		hgvs2 = pHGVS(self.ref)
		message = "%s, %s, %s, %s." % (self._testMethodName, self.query, self.ref, self.check)
		self.assertEqual(pHGVS.check_p_hgvs(self.query, self.ref ), self.check)


	#--------------------------------------------
	# check synonymous variants
	#--------------------------------------------

	def test_synonymous_1(self):
		self.query ="p.Pro559Pro"
		self.ref="p.(=)"
		self.check="yes_m"

		hgvs1 = pHGVS(self.query)
		hgvs2 = pHGVS(self.ref)
		message = "%s, %s, %s, %s." % (self._testMethodName, self.query, self.ref, self.check)
		self.assertEqual(pHGVS.check_p_hgvs(hgvs1, hgvs2 ), self.check)


	def test_synonymous_2(self):
		self.query ="p.Pro559Pro"
		self.ref="p.="
		self.check="yes_m"


		hgvs1 = pHGVS(self.query)
		hgvs2 = pHGVS(self.ref)
		message = "%s, %s, %s, %s." % (self._testMethodName, self.query, self.ref, self.check)
		self.assertEqual(pHGVS.check_p_hgvs(hgvs1, hgvs2 ), self.check)

	def test_synonymous_3(self):
		self.query ="p.Pro559="
		self.ref="p.(=)"
		self.check="yes_m"

		hgvs1 = pHGVS(self.query)
		hgvs2 = pHGVS(self.ref)
		message = "%s, %s, %s, %s." % (self._testMethodName, self.query, self.ref, self.check)
		self.assertEqual(pHGVS.check_p_hgvs(hgvs1, hgvs2 ), self.check)


	def test_synonymous_4(self):
		self.query ="p.Pro559Pro"
		self.ref="p.Pro559="
		self.check="yes_m"

		hgvs1 = pHGVS(self.query)
		hgvs2 = pHGVS(self.ref)
		message = "%s, %s, %s, %s." % (self._testMethodName, self.query, self.ref, self.check)
		self.assertEqual(pHGVS.check_p_hgvs(hgvs1, hgvs2 ), self.check)


	def test_synonymous_5(self):
		self.query ="p.Pro559="
		self.ref="p.="
		self.check="yes_m"

		hgvs1 = pHGVS(self.query)
		hgvs2 = pHGVS(self.ref)
		message = "%s, %s, %s, %s." % (self._testMethodName, self.query, self.ref, self.check)
		self.assertEqual(pHGVS.check_p_hgvs(hgvs1, hgvs2 ), self.check)


	def test_synonymous_6(self):
		self.query ="p.(=)"
		self.ref="p.Pro559="
		self.check="yes_m"

		hgvs1 = pHGVS(self.query)
		hgvs2 = pHGVS(self.ref)
		message = "%s, %s, %s, %s." % (self._testMethodName, self.query, self.ref, self.check)
		self.assertEqual(pHGVS.check_p_hgvs(hgvs1, hgvs2 ), self.check)


	def test_synonymous_7(self):
		self.query = "p.="
		self.ref="p.(=)"
		self.check="yes_m"

		hgvs1 = pHGVS(self.query)
		hgvs2 = pHGVS(self.ref)
		message = "%s, %s, %s, %s." % (self._testMethodName, self.query, self.ref, self.check)
		self.assertEqual(pHGVS.check_p_hgvs(hgvs1, hgvs2 ), self.check)


	def test_synonymous_8(self):
		self.query ="p.Pro559Glu"
		self.ref="p.(=)"
		self.check="no"


		hgvs1 = pHGVS(self.query)
		hgvs2 = pHGVS(self.ref)
		message = "%s, %s, %s, %s." % (self._testMethodName, self.query, self.ref, self.check)
		self.assertEqual(pHGVS.check_p_hgvs(hgvs1, hgvs2 ), self.check)

	def test_synonymous_9(self):
		self.query ="p.Ter257="
		self.ref="p.(=)"
		self.check="yes_m"

		hgvs1 = pHGVS(self.query)
		hgvs2 = pHGVS(self.ref)
		message = "%s, %s, %s, %s." % (self._testMethodName, self.query, self.ref, self.check)
		self.assertEqual(pHGVS.check_p_hgvs(hgvs1, hgvs2 ), self.check)


	#--------------------------------------------
	# check dups
	#--------------------------------------------

	def test_dup_1(self):
		self.query ="p.Pro559dup"
		self.ref="p.Pro559dup"
		self.check="yes"

		hgvs1 = pHGVS(self.query)
		hgvs2 = pHGVS(self.ref)
		message = "%s, %s, %s, %s." % (self._testMethodName, self.query, self.ref, self.check)
		self.assertEqual(pHGVS.check_p_hgvs(hgvs1, hgvs2 ), self.check)

	#--------------------------------------------
	# check insertions
	#--------------------------------------------

	def test_ins_1(self):
		self.query ="p.Lys2_Met3insGlnSerLys"
		self.ref="p.Lys2_Met3insGlnSerLys"
		self.check="yes"

		hgvs1 = pHGVS(self.query)
		hgvs2 = pHGVS(self.ref)
		message = "%s, %s, %s, %s." % (self._testMethodName, self.query, self.ref, self.check)
		#self.assertEqual(pHGVS.check_insertion(hgvs1, hgvs2 ), self.check)
		self.assertEqual(pHGVS.check_p_hgvs(hgvs1, hgvs2 ), self.check)



	def test_ins_2(self):
		self.query ="p.Lys2_Met3insGlnSerLysArg"
		self.ref="p.Lys2_Met3ins4"
		self.check="yes_m"

		hgvs1 = pHGVS(self.query)
		hgvs2 = pHGVS(self.ref)
		message = "%s, %s, %s, %s." % (self._testMethodName, self.query, self.ref, self.check)
		#self.assertEqual(pHGVS.check_insertion(hgvs1, hgvs2 ), self.check)
		self.assertEqual(pHGVS.check_p_hgvs(hgvs1, hgvs2 ), self.check)


	#--------------------------------------------
	# check duplication
	#--------------------------------------------

	def test_dup_p2(self):
		self.query ="p.Gly54_Arg55dupArgArg"
		self.ref="p.Gly54_Arg55dup"
		self.check="p.Gly54_Arg55dup"
		self.hgvs = pHGVS(self.query)
		self.ref = pHGVS(self.ref)

		message = "%s, %s, %s, %s." % (self._testMethodName, self.query, self.ref, self.check)
		self.assertEqual(pHGVS.check_p_hgvs(self.ref, self.hgvs.alias), 'yes_m')

	#--------------------------------------------
	# check deletions
	#--------------------------------------------

	def test_del_1(self):
		self.query ="p.Pro559del"
		self.ref="p.Pro559del"
		self.check="yes"

		message = "%s, %s, %s, %s." % (self._testMethodName, self.query, self.ref, self.check)
		self.assertEqual(pHGVS.check_p_hgvs(self.query, self.ref ), self.check)

	def test_del_2(self):
		self.query ="p.Cys28_Met30del"
		self.ref="p.Cys28_Met30del"
		self.check="yes"

		message = "%s, %s, %s, %s." % (self._testMethodName, self.query, self.ref, self.check)
		self.assertEqual(pHGVS.check_p_hgvs(self.query, self.ref ), self.check)

	def test_del_3(self):
		self.query ="p.Cys28_Met30delArgLys"
		self.ref="p.Cys28_Met30del"
		self.check="yes_m"

		message = "%s, %s, %s, %s." % (self._testMethodName, self.query, self.ref, self.check)
		self.assertEqual(pHGVS.check_p_hgvs(self.query, self.ref ), self.check)


	#--------------------------------------------
	# check insertion deletions
	#--------------------------------------------

	def test_protein_delins_1(self):
		self.query ="p.GluValThrTrp33359ValLysGluLys"
		self.ref="p.Glu33359_Trp33362delinsValLysGluLys"
		self.check="yes_m"
		message = "%s, %s, %s, %s." % (self._testMethodName, self.query, self.ref, self.check)
		self.assertEqual(pHGVS.check_p_hgvs(self.query, self.ref ), self.check)

	def test_protein_delins_1(self):
		self.query ="p.Glu22_Val23delinsLys"
		self.ref="p.Glu22_Val23delGluValinsLys"
		self.check="yes_m"
		message = "%s, %s, %s, %s." % (self._testMethodName, self.query, self.ref, self.check)
		self.assertEqual(pHGVS.check_p_hgvs(self.query, self.ref ), self.check)

	def test_protein_delins_1(self):
		self.query ="p.Gln256>ArgGlu"
		self.ref="p.Gln256delinsArgGlu"
		self.check="yes_m"
		message = "%s, %s, %s, %s." % (self._testMethodName, self.query, self.ref, self.check)
		self.assertEqual(pHGVS.check_p_hgvs(self.query, self.ref ), self.check)



	#--------------------------------------------
	# check missense variants
	#--------------------------------------------


	def test_missense_p1(self):
		self.query ="p.Pro559Ser"
		self.ref="p.Pro559Thr"
		self.check="no"

		message = "%s, %s, %s, %s." % (self._testMethodName, self.query, self.ref, self.check)
		self.assertEqual(pHGVS.check_p_hgvs(self.query, self.ref ), self.check)

	def test_missense_p2(self):
		self.query ="p.Pro559Ser"
		self.ref="p.Pro559Ser"
		self.check="yes"

		message = "%s, %s, %s, %s." % (self._testMethodName, self.query, self.ref, self.check)
		self.assertEqual(pHGVS.check_p_hgvs(self.query, self.ref ), self.check)


	#--------------------------------------------
	# check extensions
	#--------------------------------------------

	def test_check_ext_1(self):
		self.query = 'p.*463Arg'
		self.ref="p.Ter463Argext*?"
		self.check="yes_m"

		message = "%s, %s, %s, %s." % (self._testMethodName, self.query, self.ref, self.check)
		self.assertEqual(pHGVS.check_p_hgvs(self.query, self.ref ), self.check)

	def test_check_ext_2(self):
		self.query = 'p.Ter522Ser'
		self.ref="p.Ter522Serext*?"
		self.check="yes_m"

		message = "%s, %s, %s, %s." % (self._testMethodName, self.query, self.ref, self.check)
		self.assertEqual(pHGVS.check_p_hgvs(self.query, self.ref ), self.check)

	def test_check_ext_2(self):
		self.query = 'p.Ter1313Tyr'
		self.ref="p.Ter1313TyrextTer66"
		self.check="yes_m"

		message = "%s, %s, %s, %s." % (self._testMethodName, self.query, self.ref, self.check)
		self.assertEqual(pHGVS.check_p_hgvs(self.query, self.ref ), self.check)



class test_pHGVS(unittest.TestCase):

	def test_pHGVS_parse_phgvs_string(self):
		hgvs_str = 'NP_017547.3:(p.Arg34Arg)'
		var = pHGVS(hgvs_str)
		phgvs = 'p.Arg34Arg'
		refseq = 'NP_017547.3'
		print("Check %s.." % (self._testMethodName))
		self.assertEqual((phgvs, refseq), (var.name, var.refseq))

#----------------------------------------------------
# test pHGVS amino acid conversoin
#-------------------------------------------------------

# test amino acid converstion
class test_convertAminoAcidSinglet(unittest.TestCase):

	def setUp(self):
		self.longMessage=True

	def test_convertAminoAcidSinglet_1(self):
		self.hgvs_str = 'p.R652*'
		self.hgvs_triplet="p.Arg652*"
		triplet = pHGVS.hgvs_triplet_from_singlet(self.hgvs_str)
		print("Check %s.." % (self._testMethodName))
		self.assertEqual(self.hgvs_triplet, triplet )

	def test_convertAminoAcidSinglet_2(self):
		self.hgvs1="p.R642_Q643ins4"
		self.hgvs_triplet="p.Arg642_Gln643ins4"

		triplet = pHGVS.hgvs_triplet_from_singlet(self.hgvs1)
		print("Check %s.." % (self._testMethodName))
		self.assertEqual(self.hgvs_triplet, triplet )




#--------------------------------------------
# check null variants
#--------------------------------------------

class test_check_null_variant(unittest.TestCase):

	def setUp(self):
		self.longMessage=True

	def test_check_null_1(self):
		self.query = None
		self.ref="p.?"
		self.check="yes_m"

		message = "%s, %s, %s, %s." % (self._testMethodName, self.query, self.ref, self.check)
		self.assertEqual(pHGVS.check_p_hgvs(self.query, self.ref ), self.check)


	def test_check_null_2(self):
		self.query ="-"
		self.ref=""
		self.check="yes"

		message = "%s, %s, %s, %s." % (self._testMethodName, self.query, self.ref, self.check)
		self.assertEqual(pHGVS.check_p_hgvs(self.query, self.ref ), self.check)

	def test_check_null_3(self):
		self.query =""
		self.ref="p.Gly23Gly"
		self.check="no"
		message = "%s, %s, %s, %s." % (self._testMethodName, self.query, self.ref, self.check)
		self.assertEqual(pHGVS.check_p_hgvs(self.query, self.ref ), self.check)


#--------------------------------------------
# check start_lost
#--------------------------------------------

class test_check_start_lost(unittest.TestCase):

	def setUp(self):
		self.longMessage=True
		
		self.query = 'p.Met24Pro'
		self.hgvs = pHGVS(self.query)
		self.ref="p.Met1Pro"
		self.hgvs2 = pHGVS(self.ref)

		#self.check="no"
		self.check="yes"
		#self.check_hgvs = pHGVS.normalize_start_lost(self.hgvs)
		#self.check_hgvs2 = pHGVS.normalize_start_lost(self.hgvs2)
		message = "%s, %s, %s, answer: %s, %s." % (self._testMethodName, self.query, self.ref, self.check, self.hgvs.alias)

	def test_normalize_start_lost_1(self):
		self.assertEqual(self.hgvs.alias, 'p.Met24Pro')

	def test_normalize_start_lost_2(self):
		self.assertEqual(self.hgvs2.alias, 'p.Met1?')


	def test_start_lost(self):
		self.query = 'p.Met1?'
		self.ref="p.Met1Ile"
		self.check="yes_m"
		message = "%s, %s, %s, %s." % (self._testMethodName, self.query, self.ref, self.check)
		self.assertEqual(pHGVS.check_p_hgvs(self.query, self.ref ), self.check)
