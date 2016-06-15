__author__ = 'jyen'


from hgvslib import constants as c
from hgvslib.cHGVS import cHGVS
from hgvslib.pHGVS import pHGVS
import logging



def is_null(hgvs_str):
	'''
	Defines null or blank fields.
	:param hgvs_str:
	:return: True or False
	'''
	if hgvs_str in c.NULL_SET or not hgvs_str:
		return True
	else:
		return False


def check_hgvs_status(hgvs1, hgvs2):
	'''
	Primitively compares two hgvs objects for syntax. This assumes that aliases have
	been made for the hgvs objects for comparison.
	:param hgvs1: Variant object
	:param hgvs2: Variant object
	:return: yes, yes_m or no
	'''

	# for null variants, for example non-coding varaints
	if      not hgvs1 and not hgvs2:
		return c.EXACT

	# where the hgvs strings match perfectly
	elif    hgvs1.name == hgvs2.name:
		return c.EXACT

	# where normalized hgvs strings are equivalent
	elif    hgvs1.alias or hgvs2.alias:

		if (hgvs1.name == hgvs2.alias or
			hgvs1.alias == hgvs2.name or
			hgvs1.alias == hgvs2.alias):
			return c.EQUIVALENT

		else:
			return c.NO_MATCH

	# another null variant case
	elif    (is_null(hgvs1.name) and
			 is_null(hgvs2.name)):
		return c.EXACT

	# otherwise, the hgvs strings are not equivalent
	else:
		return c.NO_MATCH



def compare_hgvs(str1, str2):
	'''
	Checks whether is coding or protein hgvs from string, calls check_p_hgvs or check_c_hgvs functions
	:param str1: hgvs string
	:param str2: hgvs string
	:return: exact (yes), yes modified or equivalent (yes_m), no match (no)
	'''

	if is_null(str1) and is_null(str2):
		return c.EXACT

	if pHGVS.is_protein(str1) or pHGVS.is_protein(str2):

		hgvs1 = pHGVS(str1)
		hgvs2 = pHGVS(str2)

		return pHGVS.check_p_hgvs(hgvs1, hgvs2)

	elif cHGVS.is_coding(str1) or cHGVS.is_coding(str2):

		hgvs1 = cHGVS(str1)
		hgvs2 = cHGVS(str2)
		return cHGVS.check_c_hgvs(hgvs1, hgvs2)

	elif str1.startswith('g.') and str2.startswith('g.'):

		logging.error('Cannot compare g. hgvs %s, %s' % (str1, str2))
	else:
		logging.error('Unknown input reference:%s,%s' % (str1, str2) )
		return c.UNKNOWN

