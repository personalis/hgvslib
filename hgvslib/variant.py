__author__ = 'jyen'


from hgvslib.transcript import Transcript
from hgvslib.cHGVS import cHGVS
from hgvslib.pHGVS import pHGVS
from hgvslib.shared_functions import is_null


# Variant information for printing
VARIANT_INFO = '''
location: {loc},{ref},{alt}
ghgvs: {ghgvs}
transcript: {trans}
chgvs: {chgvs}
phgvs: {phgvs}
gene: {gene}
type: {type}
'''



class Variant(object):

	'''
	Overarching HGVS Variant Object.
	Stores multiple information about the variant, transcript, a unique identifier, cHGVS and pHGVS information.
	Creates and stores cHGVS and pHGVS objects.
	Can be used for comparing different HVS strings.
	'''

	def __init__(self, chgvs_str, phgvs_str, id, transcript='', effect='', ghgvs=''):

		# parse the chgvs string with the cHGVS object
		self.chgvs_obj  = cHGVS(chgvs_str)
		self.chgvs      = self.chgvs_obj.name
		self.type       = self.chgvs_obj.type

		self.phgvs_obj  = pHGVS(phgvs_str)
		self.phgvs      = self.phgvs_obj.name
		self.ghgvs      = ghgvs

		self.id         = id
		self.effect     = effect
		self.alt        = ''
		self.ref        = ''
		self.chr        = ''
		self.start      = ''
		self.stop       = ''
		self.location   = ''
		self.gen_alt    = ''

		self.gene       = ''

		self.transcript = transcript

		# get transcript from cHGVS object if the transcript field is null
		if not is_null(self.chgvs_obj.transcript) and is_null(transcript):
			self.transcript = self.chgvs_obj.transcript

		self.protein    = ''
		self.version    = Transcript(self.transcript).version

		self.check      = ''
		self.refalt     = ''

	@property
	def variant(self):
		return (VARIANT_INFO.format(
					loc=self.location,
					ref=self.ref,
					alt=self.alt,
					ghgvs=self.ghgvs,
					chgvs=self.chgvs,
					phgvs=self.phgvs,
					gene=self.gene,
					type=self.type
			)
		)

	def __str__(self):
		return (VARIANT_INFO.format(
					loc=self.location,
					ref=self.ref,
					alt=self.alt,
					ghgvs=self.ghgvs,
					trans=self.transcript,
					chgvs=self.chgvs,
					phgvs=self.phgvs,
					gene=self.gene,
					type=self.type
			)
		)



	@classmethod
	def from_vcf_info_field(cls, id, rec_info_field):
		'''
		Creates a Variant object from the SnpEff or VEP VCF info field (v.4.1L.
		:param id:  unique identifier to associate with teh variant
		:param rec_info_field:  VCF info field, which looks like this:
		T|stop_gained|HIGH|CXorf57|ENSG00000147231|transcript|ENST00000372548|protein_coding|4/14|c.997C>T|p.Arg333*|1106/3861|997/2568|333/855||
		:return: Variant object
		'''
		list = rec_info_field.split('|')
		alt_vcf = list[0]
		csq = list[1]
		hgvs_c = list[9]
		hgvs_p = list[10]
		transcript = list[6]

		if '=' not in hgvs_p and '(' in hgvs_p:
			hgvs_p = hgvs_p.split('(')[1].replace(')','')
			hgvs_p = hgvs_p.replace('%3D','=') # this encodes the = sign

		if hgvs_c.count(':') > 1 or hgvs_p.count(':') > 1:
			var = Variant('', '', id, '') # ignore the variant?
		else:
			var = Variant(hgvs_c, hgvs_p, id, transcript, csq)
			var.alt = alt_vcf
		return var
