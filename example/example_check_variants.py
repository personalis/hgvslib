import argparse

__author__ = 'jyen'

from hgvslib.class_functions import compare_hgvs
from hgvslib.constants import NULL_SET, NO_MATCH

def write_contents(outfile):
	try:
		return open( '%s' % outfile, 'w')
	except IOError:
		logging.error('Error: can\'t find file %s or read data' % outfile)


def check_ref_list(ref_str, hgvs_list, transcripts):
	'''
	Compares hgvs given a reference hgvs, and list of hgvs string and corresponding transcripts.
	:param ref_str: reference hgvs string
	:param hgvs_list: list of hgvs strings e.g. [c.4324A>T, c.4324A>T, c.4324A>T]
	:param transcripts: list of transcripts e.g [NM_0000342.1, NM_0000342.1, NM_0000342.1]
		*order must correspond with the hgvs strings*
	:return: list of comparison results in order e.g. [yes, yes, yes]
	'''
	r = []

	for num in [0,1,2]:

		query_str = hgvs_list[num]
		transcript = transcripts[num]

		if transcript == NO_MATCH: # only compare transcripts with exact match
			result = NO_MATCH
		else:
			result = compare_hgvs(ref_str, query_str)
		r.append(result)

	return '\t'.join(r)


def main(input_fileobject, outfile):

	# write header
	out = write_contents(outfile)
	header = input_fileobject.readline()
	column_names = header.strip('\n').split('\t')
	newheader = header.rstrip() + '\t'.join(['s_c_check','vr_c_check','vep_c_check',
	                                            's_p_check','vr_p_check','vep_p_check', '\n'])

	out.write(newheader)

	print('Reading %s, printing to %s' % (input_fileobject.name, outfile))

	for line in input_fileobject.readlines():
		line = line.rstrip().replace('NULL','')

		data = dict(zip(column_names,line.split('\t')))

		if data['transcript'] in NULL_SET:
			continue

		if all(string in NULL_SET for string in [data['snpeff_transcript'],data['vr_transcript'],data['vep_transcript']]):
			continue

		chgvs_list = [data['snpeff_c_hgvs'],data['vr_c_hgvs'],data['vep_c_hgvs']]
		phgvs_list = [data['snpeff_p_hgvs'],data['vr_p_hgvs'],data['vep_p_hgvs']]
		transcript_list = [data['snpeff_transcript'],data['vr_transcript'],data['vep_transcript']]

		chgvs_check_list = check_ref_list(data['ref_c_hgvs'], chgvs_list, transcript_list)
		phgvs_check_list = check_ref_list(data['ref_p_hgvs'], phgvs_list, transcript_list)
		
		items = '{}\t{}'.format(chgvs_check_list, phgvs_check_list)
		s = f'{line}\t{items}\n'
		out.write(s)

	out.close()


if __name__ == '__main__':

	parser = argparse.ArgumentParser(description='Parse compare output for hgvs checking ')
	parser.add_argument('-i', '--infile', nargs='?', type=argparse.FileType('r'), required = True)
	args = parser.parse_args()
	outfile = str(args.infile.name)[:-3] + 'parsed.txt'
	main(args.infile, outfile)

