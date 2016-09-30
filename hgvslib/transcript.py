__author__ = 'jyen'

from hgvslib import constants as c

class Transcript(object):
	'''
	Transcript class. Used for comparing transcript versions and accessions.
	'''
	def __init__(self, name):

		self.name = name
		self.version = ''
		if '.' in self.name:
			self.accession, self.version = self.name.split('.')
		else:
			self.version = ''
			self.accession = self.name

	# checks if the transcript matches a given transcript input
	def check_transcript_ref(self,ref):

		print "checking", self.name, ref
		if self.name == ref:
			print "exact"
			return c.EXACT

		elif self.name in c.NULL_SET or ref in c.NULL_SET or not ref:
			print "null"
			return c.NO_MATCH

		elif '\.' in ref:
			ref, ref_version = ref.split('\.')
			print "looking at version",self.accession, ref
			if self.accession == ref:
				print "version"
				return 'yes,version'
			else:
				print "no match version"
				return c.NO_MATCH

		else:
			print "no match"
			return c.NO_MATCH
