__author__ = 'jyen'

import re
import logging

from hgvslib import constants as c
from hgvslib.shared_functions import is_null, check_hgvs_status

from hgvslib.pythonalis.aminoacid import AminoAcid

#---------------------------------------------------------------
# PROTEIN HGVS
#---------------------------------------------------------------

class pHGVS(object):

	'''
	pHGVS object.
	Normalized variant expression, get variant type.
	Used for comparing pHGVS nomenclature.
	'''

	def __init__( self, hgvs_str ):
		refseq_id, hgvs_str = pHGVS.parse_phgvs_string(hgvs_str)
		self.refseq     = refseq_id
		self.name       = hgvs_str
		self.alias      = hgvs_str   # default
		self.pos        = ''  # position of change
		self.aa1        = ''  # first amino acid
		self.aa2        = ''  # second amino acid

		self._get_pHGVS_type()
		self._get_normalized_alias()


	def _get_pHGVS_type(self):
		'''
		Return the variant type (e.g. effect impact) based on variant syntax.
		Note: order of the statements is important here. For example, p.Arg1921MetfsTer9* is
		foremost a frameshift variant, than nonsense variant.
		:return: variant type based on pHGVS syntax
		'''


		if self.name in c.NULL_SET:
			self.type = '?'

		elif 'p.' not in self.name and self.name not in c.NULL_SET:
			raise Exception('Not a protein coding variant.')

		elif self.__is_synonymous():
			self.type = c.SYNONYMOUS

		elif self.__is_indel():
			self.type = c.DELINS

		elif c.DEL in self.name and c.INS not in self.name:
			self.type = c.DEL

		elif c.DUP in self.name:
			self.type = c.DUP

		elif c.FRAMESHIFT in self.name:
			self.type = c.FRAMESHIFT

		elif self.__is_nonsense():
			self.type = c.NONSENSE

		elif c.INS in self.name:
			self.type = c.INS

		elif self.__is_start_lost():
			self.type = c.START_LOST

		elif self.name in c.NULL_SET or not self.name:
			self.type = c.NON_CODING

		elif self.__is_extension():
			self.type = c.EXTENSION

		elif '_' in self.name and '>' in self.name:
			self.type = c.DINUCLEOTIDE

		else:
			self._get_snv_type()


	def __is_indel(self):
		'''
		Checks if variant is an indel.
		'''

		if (c.DELINS in self.name or '>' in self.name or
			(c.DEL in self.name and c.INS in self.name)):
			# frameshift variant takes precedent
			if 'fs>' not in self.name and 'fs*>' not in self.name:
				return True
		else:
			return False

	def __is_start_lost(self):
		'''
		Checks if variant is in the start codon.
		'''

		if (self.name.startswith('p.Met1') or
			  self.name == 'p.Met?'):
			return True
		else:
			return False

	def __is_extension(self):
		'''
		Checks if variant is a terminator codon.
		'''
		hgvs = self.name.replace('Ter','*')
		if c.EXTENSION in hgvs or \
				(hgvs.startswith('p.*') and not
				hgvs.endswith('*')):
			return True
		else:
			return False

	def __is_nonsense(self):
		'''
		Checks if variant is nonsense.
		:return: True/False
		'''
		hgvs = self.name.replace('Ter','*')
		if ( not hgvs.startswith('p.*') and hgvs.endswith('*') ):
			return True
		else:
			return False

	def __is_synonymous(self):
		'''
		Checks if variant is synonymous. If it's synonymous also assigns the alias to p.(=)
		:return: True/False
		'''

		if (self.name == c.HGVS_SYN or
			self.name in c.SYN_ALIAS_SET or
			self.name.endswith('=')):
			self.alias = c.HGVS_SYN
			return True
		else:
			# need to get the amino acids for synonymous cases that are p.Gly235=
			hgvs_re = re.search(r'p[.]([a-zA-Z]+)(\d+)([a-zA-Z]+)', self.name)
			if hgvs_re:
				aa1, pos, aa2 = hgvs_re.group(1, 2, 3)
				if aa1 == aa2:
					self.alias = c.HGVS_SYN
					self.aa1 = aa1
					self.aa2 = aa2
					self.pos = pos
					return True
				else:
					return False
			else:
				return False

	def _get_snv_type(self):
		'''
		Checks if variant is an SNV or multi SNV.
		 e.g. p.Thr29Thr vs p.GluValThrTrp33359ValLysGluLys
		Sets the amino acids to Variant object.
		:return: missense, multi_sub or unknown.
		'''

		hgvs_re = re.search(r'p[.]([a-zA-Z]+)(\d+)([a-zA-Z]+)', self.name)
		hgvs_re2 = re.search(r'p[.]([*])(\d+)([*])', self.name.replace('Ter',('*')))

		if hgvs_re:
			aa1, pos, aa2 = hgvs_re.group(1, 2, 3)
			if aa1 == aa2:
				self.type = c.SYNONYMOUS

			# for p.Thr29Pro
			elif len(aa1) == 3 and len(aa2) == 3:
				self.aa1 = aa1
				self.aa2 = aa2
				self.pos = pos
				self.type = c.MISSENSE

			# for p.GluValThrTrp33359ValLysGluLys
			elif len(aa1)/3 == len(aa2)/3:
				self.aa1 = aa1
				self.aa2 = aa2
				self.pos = pos
				self.type = c.MULTI_SUB
			else:
				self.type = c.UNKNOWN

		elif hgvs_re2:
			aa1, pos, aa2 = hgvs_re2.group(1, 2, 3)
			self.type = c.SYNONYMOUS
			self.pos = pos
			self.alias = c.HGVS_SYN

		else:
			self.type = c.UNKNOWN

	def _get_normalized_alias(self):
		'''
		Normalize the syntax for the pHGVS variants to its most minimal form.
		:return: a normalized, minimal syntax
		'''

		if self.type is c.INS:
			self._normalize_insertion()

		elif self.type is c.EXTENSION:
			self._normalize_extension()

		elif self.type is c.NONSENSE:
			self._normalize_nonsense()

		elif self.type is c.FRAMESHIFT:
			self._normalize_frameshift()

		elif self.type is c.DEL:
			self._normalize_del()

		elif self.type is c.DELINS:
			self._normalize_delins()

		elif self.type is c.START_LOST:
			self._normalize_start_lost()

		elif self.type is c.MULTI_SUB:
			self._normalize_multi_sub()

	#---------------------------------------------------------------
	# The following are functions to normalize pHGVS syntax
	#---------------------------------------------------------------
	def _check_amino_acids(self):
		'''
		Stores the amino acids and position from the variant syntax.
		E.g. p.Arg222= will store aa1 and aa2 as Arg, and pos as 222.
		Must not be  synonymous: e.g. p.(=) or p.= - . Need to have amino acids for this function
		:return: nothing
		'''
		if self.name.endswith('='): # must not be of p.(=) type
			hgvs_re  = re.search(r'p[.]([a-zA-Z]+)(\d+)[=]', self.name)
			if hgvs_re:
				self.aa1, self.pos = hgvs_re.group(1, 2)
				self.aa2 = self.aa1

		else:  # this is for p.Thr29Thr
			hgvs_re = re.search(r'p[.]([a-zA-Z]+)(\d+)([a-zA-Z]+)', self.name)

			if hgvs_re:
				self.aa1, self.pos, self.aa2 = hgvs_re.group(1, 2, 3)

	def _normalize_insertion(self):
		'''
		Normalize insertion to its most minimal form.
			e.g. p.Arg54_Gly55insGluArgGlu to p.Arg54_Gly55ins3
		:return: Syntax so that the number of amino acids are represented at the end.
		'''
		aa_list = self.name.split(c.INS)

		if not aa_list[1].isdigit():
			num_aa = len(aa_list[1])/3
			self.alias = '{}{}{}'.format(aa_list[0], c.INS, num_aa)

	def _normalize_extension(self):
		'''
		Normalize extension p.*1258Tyr format or p.Ter1258Tyrext*?  to p.*1258Tyrext*?
		p.Ter346Serext*?,
		'''
		self.alias = self.name.replace('Ter','*')

		# check if end has extension
		parts = self.alias.split(c.EXTENSION)

		last_aa = self.name[-3:]
		if last_aa in c.AMINO_ACID_TRIPLETS:
			self.alias += 'ext*?'
		else:
			self.alias = parts[0] + 'ext*?'

	def _normalize_nonsense(self):
		'''
		Normalize termination to * instead of Ter
		e.g. p.Arg222Ter to p.Arg222*
		:return:
		'''
		self.alias = self.name.replace('Ter','*')

	def _normalize_frameshift(self):
		'''
		Convert frameshift variants to short form p.hgvs
		e.g. p.Glu67fs*10 or p.Glu67Glyfs* to p.Glu67fs
		:return:
		'''
		if c.FRAMESHIFT in self.name:
			self.alias = self.name.split(c.FRAMESHIFT)[0] + c.FRAMESHIFT

		# convert p.Glu67Glyfs* back to p.Glu67fs
		hgvs_re = re.search(r'p[.]([A-Z][a-z][a-z])(\d+)([A-Z][a-z][a-z])fs', self.name)

		if hgvs_re:
			first_amino, pos, second_amino = hgvs_re.group(1, 2, 3)
			self.alias = 'p.{}{}fs'.format(first_amino, pos)

	def _normalize_del(self):
		'''
		Normalize deletion syntax
		e.g. p.Arg222del vs p.Arg222del1 vs p.Arg222delA

		:return: normalized deletion syntax
		'''
		if self.type is c.DEL:
			self.alias = self.name.split(c.DEL)[0] + c.DEL

	def _normalize_start_lost(self):
		'''
		Normalize start lost syntax:
		e.g. p.Met1Lys to p.Met1?
		'''
		if self.type is c.START_LOST:
			hgvs_re = re.search(r'p[.]Met1([A-Z][a-z][a-z])', self.name)
			if hgvs_re:
				amino = hgvs_re.group(1)
				if amino in c.AMINO_ACID_TRIPLETS:
					self.alias= 'p.Met1?'

	def _normalize_duplication(self):
		'''
		Normalize duplication syntax.
		This function can only be used if there is evidence that it may be a dup. Because we are not actually checking
		the genomic reference that the aas at this position match the inserted aas
		e.g. p.Arg222_Arg223insArgArg  to p.Arg222_Arg223dup
		:return:
		'''
		aa_list = self.name.split(self.type)
		if aa_list:
			#  p.Arg222_Arg223dup vs  p.Arg222_Arg223insArgArg -> if the second part is not a digit, convert to number of amino acids
			hgvs_re = re.search(r'p[.]([A-Z][a-z][a-z])(\d+)_([A-Z][a-z][a-z])(\d+)', aa_list[0])
			if hgvs_re:
				first_amino, first_pos, second_amino, second_pos = hgvs_re.group(1, 2, 3, 4)
				if first_amino == second_amino:
					if not aa_list[1].isdigit():
						aas = aa_list[1] # the amino acid string
						len_aa = len(aas)/3 # the number of amino acids
						# if the inserted amino acids are the same e.g. p.Arg54_Arg55insArgArg
						if aas == aas[:3]*len_aa and first_amino == aas[:3]:
							self.alias = self.name.split(self.type)[0] + c.DUP

	def _normalize_delins(self):
		'''
		Reformats indels to their minimal forms
		e.g. p.Arg222_Glu223delArgGluinsGly to p.Arg222_Glu223delinsGly
		e.g. for COSMIC indel annotation - for p.Glu23>GlyArg, return p.Glu23delinsGlyArg
		'''
		if '>' in self.name:
			self.alias = self.name.replace('>',c.DELINS)
		else:
			try:
				start = self.name.split(c.DEL)[0]
				end = self.name.split(c.INS)[1]
				self.alias = start + c.DELINS + end
			except:
				self.alias = self.name

	def _normalize_multi_sub(self):
		'''
		#Normalize start lost syntax:
		#e.g. p.Met1Lys to p.Met1?
		#e.g. p.GluValThrTrp33359ValLysGluLys vs p.Glu33359_Trp33362delinsValLysGluLys
		'''
		num_aa = len(self.aa1)/3
		self.alias = 'p.{first_aa}{first_pos}_{end_aa}{end_pos}delins{aas}'.format(
			first_aa = self.aa1[0:3],
			first_pos = self.pos,
			end_aa = self.aa1[-3:],
			end_pos = int(self.pos) + ( num_aa - 1),
			aas = self.aa2)

	@classmethod
	def parse_phgvs_string(cls, hgvs_str):
		"""
		Cleans up the hgvs string, parses it into refseq id and phgvs
		:param hgvs_str: input hgvs string
		:return: refseq accession, phgvs
		"""

		if hgvs_str is None or not hgvs_str:
			return '', ''

		elif c.PROTEIN_START not in hgvs_str and hgvs_str not in c.NULL_SET:
			logging.error('Not a protein coding variant: %s.' % hgvs_str)
			return '', ''

		# some formatting
		elif ':' in hgvs_str:
			refseq_id, phgvs = hgvs_str.split(':')
			# remove brackets that are not associated with p.(=)
			if phgvs.startswith('(') and phgvs.endswith(')'):
				phgvs = phgvs.rstrip(')').lstrip('(')
			return refseq_id, phgvs

		else:
			# remove brackets that are not associated with p.(=)
			if hgvs_str.startswith('(') and hgvs_str.endswith(')'):
				hgvs_str = hgvs_str.rstrip(')').lstrip('(')
			return '', hgvs_str


	#---------------------------------------------------------------
	# Comparing pHGVS syntax
	#---------------------------------------------------------------

	@classmethod
	def is_phgvs_instance(cls, hgvs_obj1):
		'''
		Creates a new instance of hgvs_obj1 if input is string. Hgvs_obj2 is only for comparison
		:param hgvs_obj1:
		:param hgvs_obj2:
		:return:
		'''
		if not isinstance(hgvs_obj1, pHGVS):
			if type(hgvs_obj1) is str:
				if is_null(hgvs_obj1):
					return ''
				else:
					return pHGVS(hgvs_obj1)
			else:
				logging.error('Not pHGVS object %s' % hgvs_obj1 )
				return hgvs_obj1
		else:
			return hgvs_obj1


	@classmethod
	def check_p_hgvs(cls, hgvs_obj1, hgvs_obj2):


		# ensure that strings are converted into pHGVS instances except for null variants
		hgvs_obj1 = pHGVS.is_phgvs_instance(hgvs_obj1)
		hgvs_obj2 = pHGVS.is_phgvs_instance(hgvs_obj2)

		# check to see if both input hgvs_obj are null variants
		if (is_null(hgvs_obj1) and is_null(hgvs_obj2)):
			if hgvs_obj1 == hgvs_obj2:      return c.EXACT
			else:                           return c.EQUIVALENT

		# basically because an empty hgvs_obj will not have a name attribute, need to make sure that the
		# pHGVS object exists for both hgvs_obj1 and hgvs_obj3, and then check that the name strings exist
		elif (not (not hgvs_obj1) and not (not hgvs_obj2)) and (not hgvs_obj1.name and not hgvs_obj2.name):
			if hgvs_obj1 == hgvs_obj2:      return c.EXACT
			else:                           return c.EQUIVALENT			

		elif is_null(hgvs_obj1) or is_null(hgvs_obj2):
			return c.NO_MATCH

		elif hgvs_obj1.type is c.DUP or hgvs_obj2.type is c.DUP:
			return pHGVS.check_dup(hgvs_obj1, hgvs_obj2)

		else:
			result = check_hgvs_status(hgvs_obj1, hgvs_obj2)
			if not result:
				print 'ERROR in result', hgvs_obj1.name, hgvs_obj2.name
			return result


	@classmethod
	def check_dup(cls, hgvs1_obj, hgvs2_obj):
		'''
		If the combination is an insertion variant and duplication, check if insertion is really a duplication.
		e.g. p.Arg222_Arg223dup vs  p.Arg222_Arg223insArgArg
		Otherwise convert to their minimal forms and check if duplications are equivalent.

		:param hgvs1_obj: pHGVS object
		:param hgvs2_obj: pHGVS object
		:return: exact, equivalent, or no match
		'''

		if hgvs1_obj.type is c.INS and hgvs2_obj.type is c.DUP:
			hgvs1_obj._normalize_duplication()

		elif hgvs2_obj.type is c.INS and hgvs1_obj.type is c.DUP:
			hgvs2_obj._normalize_duplication()

		else:
			hgvs1_obj.alias = hgvs1_obj.name.split(hgvs1_obj.type)[0] + hgvs1_obj.type
			hgvs2_obj.alias = hgvs2_obj.name.split(hgvs2_obj.type)[0] + hgvs2_obj.type

		return check_hgvs_status(hgvs1_obj, hgvs2_obj)



	@classmethod
	def replace_amino_acid_singlet(cls, amino_singlet_re):
		'''
		Finds all single amino acid occurrences and replaces them with triplet amino acid
		'''
		amino_singlet = None
		if amino_singlet_re is not None:
			amino_singlet = amino_singlet_re.group(1)
		if amino_singlet and amino_singlet in c.AMINO_ACID_SINGLETS:
			return AminoAcid.triplet_from_singlet(amino_singlet)
		return amino_singlet

	@classmethod
	def replace_amino_acid_triplet(cls, amino_acid_triplet_re):
		'''
		Finds all triplet amino acid occurrences and replaces them with singlet amino acid
		'''
		amino_acid_triplet = None
		if amino_acid_triplet_re is not None:
			amino_acid_triplet = amino_acid_triplet_re.group(1)
		if amino_acid_triplet and amino_acid_triplet in c.AMINO_ACID_TRIPLETS:
			return AminoAcid.singlet_from_triplet(amino_acid_triplet)
		return amino_acid_triplet


	@classmethod
	def hgvs_triplet_from_singlet(cls, hgvs_str):
		'''
		Converts singlet amino acid to triplet
		'''
		if not hgvs_str:
			return ''
		hgvs = re.sub(r'([A-Z])', pHGVS.replace_amino_acid_singlet, hgvs_str)
		return hgvs

	@classmethod
	def hgvs_singlet_from_triplet(cls, hgvs_str):
		'''
		Converts singlet amino acid to triplet
		'''
		if not hgvs_str:
			return ''
		hgvs = re.sub(r'([A-Z][a-z][a-z])', pHGVS.replace_amino_acid_triplet, hgvs_str)
		return hgvs

	@classmethod
	def is_protein(cls, hgvs_str):
		'''
		Defines protein variants.
		:param hgvs_str: hgvs expression
		:return: True or False
		'''
		if (hgvs_str.startswith(c.PROTEIN_START) or
		   c.PROTEIN_START in hgvs_str):
			return True
		else:
			return False
