import argparse

__author__ = 'jyen'

import sys, os
sys.path.append(os.path.join(os.path.dirname(__file__), '..'))
from hgvslib.class_functions import compare_hgvs
from hgvslib.constants import NULL_SET, NO_MATCH
from hgvslib.pHGVS import pHGVS
from hgvslib.transcript import Transcript

def write_contents(outfile):
	try:
		return open( '%s' % outfile, 'w')
	except IOError:
		logging.error('Error: can\'t find file %s or read data' % outfile)


def check_transcript(ref_transcript_str, transcript_list):
	'''
	Compares an input reference transcript with each query transcript in a list
	:param ref_transcript_str: reference transcript
	:param transcript_list: list of transcripts to compare
	:return: list of results in order e.g. [yes, yes, version]l exact match (yes), different version (yes,version) or no match (no).
	'''
	r = []
	for num in range(0,len(transcript_list)):

		transcript = transcript_list[num]

		trans_obj = Transcript(transcript)
		check_transcript = trans_obj.check_transcript_ref(ref_transcript_str)
		#print 'transcript:', transcript, 'trans obj accession:', trans_obj.accession, "ref:", ref_transcript, check_transcript
		r.append(check_transcript)

	return '\t'.join(r)


def check_ref_list(ref_str, hgvs_list):
	'''
	Compares hgvs given a reference hgvs, and list of hgvs string and corresponding transcripts.
	:param ref_str: reference hgvs string
	:param hgvs_list: list of hgvs strings e.g. [c.4324A>T, c.4324A>T, c.4324A>T]
	:return: list of comparison results in order e.g. [yes, yes, yes]
	'''
	r = []

	for num in range(0,len(hgvs_list)):

		query_str = hgvs_list[num]
		result = compare_hgvs(ref_str, query_str)
		r.append(result)

	return '\t'.join(r)


def main(input_fileobject, outfile):

	# write header
	out = write_contents(outfile)
	header = input_fileobject.readline()
	column_names = header.strip('\n').split('\t')
	newheader = header.rstrip() + '\t' + '\t'.join(['transcript_check', 'c_check','p_check','\n'])
	out.write(newheader)

	print 'Reading %s, printing to %s' % (input_fileobject.name, outfile)

	for line in input_fileobject.readlines():
		line = line.rstrip().replace('NULL','')

		data = dict(zip(column_names,line.split('\t')))
		# skip variants with no transcript identifiers
		if data['ref_transcript'] in NULL_SET or 'test_transcript' not in data:
			continue

		if all(string in NULL_SET for string in [data['test_transcript']]):
			continue

		chgvs_list = [data['test_c_hgvs']]
		transcript_list = [data['test_transcript']]
		
		# if phgvs is in singlet form, convert to triplet for comparison
		phgvs = pHGVS.hgvs_triplet_from_singlet("p." + data['test_p_hgvs'])
		phgvs_list = [phgvs]


		chgvs_check_list = check_ref_list(data['ref_c_hgvs'], chgvs_list)
		phgvs_check_list = check_ref_list(data['ref_p_hgvs'], phgvs_list)
		transcript_check_list = check_transcript(data['ref_transcript'], transcript_list)
		items = '{}\t{}\t{}'.format(transcript_check_list, chgvs_check_list, phgvs_check_list)
		s = '{line}\t{hgvs_items}\n'.format(line=line, hgvs_items=items)
		out.write(s)

	out.close()
	print 'Done.'

if __name__ == '__main__':

	parser = argparse.ArgumentParser(description='Parse compare output for hgvs checking ')
	parser.add_argument('-i', '--infile', nargs='?', type=argparse.FileType('r'), required = True)
	args = parser.parse_args()
	outfile = str(args.infile.name)[:-3] + 'parsed.txt'
	main(args.infile, outfile)

