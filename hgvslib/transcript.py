__author__ = 'jyen'

from hgvslib import constants as c

class Transcript(object):
	'''
	Transcript class. Used for comparing transcript versions and accessions.
	'''
	def __init__(self, name):

		self.name = name
		self.no_version = ''
		if '.' in self.name:
			self.no_version, self.version = self.name.split('.')
		else:
			self.version = ''

	# checks if the transcript matches a given transcript input
	def check_transcript_ref(self,ref):

		if self.name == ref:
			return c.EXACT

		elif self.name in c.NULL_SET or ref in c.NULL_SET or not ref:
			return c.NO_MATCH

		elif '.' in ref:
			ref, ref_version = ref.split('.')
			if self.refseq == ref:
				return 'yes,version'
			else:
				return c.NO_MATCH

		else:
			return c.NO_MATCH
