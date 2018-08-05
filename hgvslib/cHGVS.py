__author__ = 'jyen'

from hgvslib import constants as c
from hgvslib.shared_functions import check_hgvs_status


import logging


# to print the variant information
VARIANT_INFO = '''
transcript: {trans}
name: {chgvs}
alias: {alias}
type: {type}
is intronic: {intronic}

'''


class cHGVS(object):
	'''
	cHGVS object.
	input: coding hgvs string in the form 'NM_00113323.1:c.232A>T'
	output: a cHGVS object with transcript and name (c.232A>T).
	Can add other attributes to the object after initialization.
	'''

	def __init__( self, hgvs_str):

		refseq, chgvs = cHGVS.parse_chgvs_string(hgvs_str)
		self.name = chgvs
		self.transcript = refseq
		self.bases          = ''
		self.type           = ''
		self.alias          = ''
		self.intronic       = ''
		self._is_intronic()
		self._get_chgvs_type()
		self._normalize_chgvs()

	def __str__(self):

		return (VARIANT_INFO.format(
			trans=self.transcript,
			chgvs=self.name,
			alias=self.alias,
			type=self.type,
			intronic=self.intronic
			)
		)



	# boolean result for whether the variant is in intronic position
	def _is_intronic(self):
		if '+' in self.name or '-' in self.name:
			self.intronic = True
		else:
			self.intronic = False

	# NM_000602.4:c.-820_-817G(4_5)
	def _get_chgvs_type(self):

		if self.name in c.NULL_SET:
			self.type =  c.UNKNOWN

		elif (c.DEL in self.name) and (c.INS in self.name):
			if c.DELINS not in self.name:
				self.type = c.DELINS_SPLIT
			else:
				self.type = c.DELINS
			self.type =  c.DELINS

		elif '_' in self.name and '>' in self.name:
			self.type =  c.DELINS

		elif c.DEL in self.name:
			self.type =  c.DEL

		elif c.INS in self.name:
			self.type =  c.INS

		elif c.DUP in self.name:
			self.type =  c.DUP

		elif '>' in self.name:
			self.type =  c.SUB

		elif self.name.endswith('='):
			self.type =  c.SYNONYMOUS

		elif self.name.startswith('*'):
			self.type =  c.UPSTREAM

		elif c.INVERSION in self.name:
			self.type =  c.INVERSION

		# for these NM_001007026.1:c.1462_1464CAG(6_35)
		elif '_' in self.name:
			self.type = c.DEL

		else:
			self.type =  c.UNKNOWN


	def _normalize_chgvs(self):
		# normalizes chgvs abbreviations to their most minimalist forms

		if self.type is c.DELINS:
			self._reformat_delins()

		elif self.type is c.DEL or self.type is c.DUP:
			self._reformat_del_or_dup()

		elif self.type is c.INS:
			self._reformat_insertion()

	def _reformat_del_or_dup(self):
		'''
		reformat deletions and duplications in alias
		stores bases for comparison later

		e.g. c.22_23dupA has alias c.22_23dup; stores 'A'
		e.g. c.22_23delAA has alias c.22_23del; stores 'AA'
		'''
		hgvs_list = self.name.split(self.type)
		self.alias = hgvs_list[0] + self.type
		try:
			self.bases = hgvs_list[1]
		except IndexError:
			pass

	def _reformat_delins(self):
		'''
		Reformats indels to their minimal forms
		e.g. convert c.222_223delTTinsGA to c.222_223delinsGA
		e.g. for COSMIC indel annotation - for c.222_223TT>GA, return c.22_223delinsGA
		e.g. for c.222_223delTTinsGA, return c.22_223delinsGA
		'''

		if '_' in self.name and '>' in self.name:
			position = self.name.split('>')[0].replace('A','').replace('G','').replace('C','').replace('T','')
			self.alias = position + c.DELINS + self.name.split('>')[1]
		else:
			self.alias = self.name.split(c.DEL)[0] + c.DELINS + self.name.split(c.INS)[1]

	def _reformat_insertion(self):
		'''
		reformats insertions to their minimal forms
		# e.g. convert c.764_765insAACCTGACAGTTGCAGTTTTCACCCATGGAAAG to c.764_765ins33
		'''
		try:
			ins_bases = self.name.split(c.INS)[1]
			if not ins_bases.isdigit():
				self.alias = '{}ins{}'.format(self.name.split(c.INS)[0], len(ins_bases))
		except IndexError:
			logging.error('Cannot normalize insertion %s' % self.name)



	@classmethod
	def parse_chgvs_string(cls, hgvs_str):
		'''
		Standalone function for parsing HGVS strings.
		:param hgvs_str: an hgvs string
		:return: separated strings for 1) accession and 2) hgvs nomenclature
		'''

		refseq = '-'
		hgvs = '-'
		# for this use case NM_001530.3:c.222_223delinsGA (c.223G>A)
		if ' ' in hgvs_str:
			hgvs_str = hgvs_str.split(' ')[0]


		# for use with ClinVar output:  NM_017547.3(FOXRED1):c.694C>T
		if c.CODING_START in hgvs_str and '(' in hgvs_str and '):c.' in hgvs_str:
			list = hgvs_str.split('(')
			refseq = list[0]
			try:
				hgvs = list[1].split(':')[1]
			except IndexError:
				logging.error('Cannot split this string %s (from %s)' % (list[1], hgvs_str) )
				pass

		# for variants with a refseq accession
		elif ':' in hgvs_str:
			try:
				list = hgvs_str.split(':')
				refseq = list[0]
				hgvs = list[1]
			except ValueError:
				logging.error('Error in hgvs string:', hgvs_str)
		# for variants with no refseq accession e.g. c.123A>T
		else:
			hgvs = hgvs_str

		return refseq, hgvs



	#----------------------------------------------------------------
	# Compare coding hgvs functions
	#----------------------------------------------------------------
	@classmethod
	def check_c_hgvs(cls, hgvs_obj1, hgvs_obj2):
		'''
		:param hgvs_obj1:  hgvs query 1
		:param hgvs_obj2: hgvs query 2
		:return: whether the two hgvs syntaxes are the same, equivalent or do not match
		'''
		# for cases where the input is not a chgvs object, create a new cHGVS object
		hgvs_obj1 = cHGVS.check_chgvs_instance(hgvs_obj1)
		hgvs_obj2 = cHGVS.check_chgvs_instance(hgvs_obj2)

		# if both are nulls
		if not hgvs_obj1 and not hgvs_obj2:
			return c.CANNOT_ASSESS

		elif hgvs_obj1.type is c.UNKNOWN or hgvs_obj2.type is c.UNKNOWN:
			return c.NO_MATCH

		# deletions and dups can be matched on position only
		elif hgvs_obj1.type is c.DUP or hgvs_obj1.type is c.DEL:
			return cHGVS.check_c_del_dup(hgvs_obj1, hgvs_obj2)

		else:
			return check_hgvs_status(hgvs_obj1, hgvs_obj2)



	@classmethod
	def check_c_del_dup(cls, hgvs1, hgvs2):
		'''
		:param hgvs1: c_hgvs deletion or duplication
		:param hgvs2: c_hgvs deletion or duplication
		:return: whether the two syntaxes are the same, equivalent or do not match
		'''
		# compares deletions or duplications
		# check if number of bases correspond e.g. c.23_25del2 vs c.23_25delCA
		# re-add bases to the alias
		try:
			if hgvs1.bases and hgvs2.bases:
				if hgvs1.bases.isdigit() and not hgvs2.bases.isdigit():
					hgvs2.alias += str(len(hgvs2.bases))
					hgvs1.alias += hgvs1.bases
				elif not hgvs1.bases.isdigit() and hgvs2.bases.isdigit():
					hgvs1.alias += str(len(hgvs1.bases))
					hgvs2.alias += hgvs2.bases
				else:
					hgvs1.alias += hgvs1.bases
					hgvs2.alias += hgvs2.bases

		except IndexError:
			pass

		return check_hgvs_status(hgvs1, hgvs2)

	@classmethod
	def check_chgvs_instance(cls, input_hgvs):
		if not isinstance(input_hgvs, cHGVS):
			if type(input_hgvs) is str:
				if (input_hgvs in c.NULL_SET or not input_hgvs):
					return ''
				else:
					return cHGVS(input_hgvs)
			else:
				logging.error('Not cHGVS object %s' % input_hgvs )
		else:
			return input_hgvs

	@classmethod
	def is_coding(cls, hgvs_str):
		'''
		Defines coding variants as transcript variants that are mitochondrial, ribosomal, nucleotide.
		Technically these might not all be 'coding' but are transcript based.
		:param hgvs_str: hgvs expression
		:return: True or False
		'''
		if (c.PROTEIN_START not in hgvs_str and
			c.GENOMIC_START not in hgvs_str and
			any(substring in hgvs_str for substring in c.CODING_START_LIST)):
			#any(substring in hgvs_str for substring in ['c.','n.','m.','r.'])):
			return True
		else:
			return False


	def get_alias(chgvs_str):
        	"""
		Converts cHGVS to the minimal form. 
		:param hgvs_str: protein HGVS format
		:return:  pHGVS in most minimal form
        	"""
        	try:
                	return(cHGVS(chgvs_str).alias)
        	except:
                	return('none')
