__author__ = 'jyen'

from hgvslib import constants as c

class Effect(object):
	'''
	Effect object.
	input: effect string (e.g. intron_variant or missense_variant)
	output: an Effect object with the effect string as the name.
	'''

	def __init__( self, name, var_type ):
		self.name       = name

	@classmethod
	def check_multi_effect(cls, effect_str, ref_str):
		'''
		To check the effect in multi-effect strings typical of VEP.
		:param effect_str: multi-effect string, e.g. intron_variant&splice_acceptor_variant
		:param ref_str: reference string, e.g. splice_acceptor_variant
		:return: yes, yes_m or no
		'''
		final_check = c.NO_MATCH
		effect_str = effect_str.replace('+','&').replace(',','&')
		for eff in effect_str.split('&'):
			final_check = Effect.check_effect(eff, ref_str)
			if final_check in [c.EXACT, c.EQUIVALENT]:
				break
		return final_check


	@classmethod
	def check_effect(cls, effect_str, ref_str):
		'''
		Checks the equivalence of effect strings.
		:param effect_str: e.g. splice_acceptor_varian
		:param ref_str: reference string, e.g. splice_acceptor_variant
		:return: yes, yes_m or no
		'''
		effect_str = effect_str.rstrip(' ').replace(' ', '_')
		ref_str = ref_str.rstrip(' ')

		if effect_str == ref_str:
			return c.EXACT
		elif effect_str in ref_str:
			return c.EQUIVALENT
		elif 'frameshift' in ref_str and 'frameshift' in effect_str:
			return c.EQUIVALENT
		else:
			for set in c.EFFECT_EQ_LIST:
				if ref_str in set and effect_str in set:
					return c.EQUIVALENT
			return c.NO_MATCH


	@classmethod
	def get_effect(cls, effect_str):
		'''
		Returns normalized effect string names.
		:param effect_str: e.g. non_synonymous_codon
		:return: normalized term e.g. missense_variant
		'''
		dict = {
			'non_synonymous_codon' : 'missense_variant',
			'upstream_variant' : 'upstream_gene_variant',
			'synonymous_codon' : 'synonymous_variant',
			'disruptive_inframe_deletion' : 'inframe_deletion',
		}

		if effect_str in dict:
			return dict[effect_str]
		else:
			return effect_str