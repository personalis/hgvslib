__author__ = 'jyen'


from hgvslib import constants as c

def is_null(hgvs_str):
	"""
	Defines null or blank fields.
	:param hgvs_str:
	:return: True or False
	"""
	if hgvs_str in c.NULL_SET or not hgvs_str:
		return True
	else:
		return False



def check_hgvs_status(hgvs1, hgvs2):
	"""
	Primitively compares two hgvs objects for syntax. This assumes that aliases have
	been made for the hgvs objects for comparison.
	:param hgvs1: Variant object
	:param hgvs2: Variant object
	:return: yes, yes_m or no
	"""

	if      not hgvs1 and not hgvs2:
		return c.EXACT

	elif    hgvs1.name == hgvs2.name:
		return c.EXACT

	elif    hgvs1.alias or hgvs2.alias:

		if (hgvs1.name == hgvs2.alias or
			hgvs1.alias == hgvs2.name or
			hgvs1.alias == hgvs2.alias):
			return c.EQUIVALENT

		else:
			return c.NO_MATCH

	elif    (is_null(hgvs1.name) and
			 is_null(hgvs2.name)):
		return c.EXACT

	else:
		return c.NO_MATCH

