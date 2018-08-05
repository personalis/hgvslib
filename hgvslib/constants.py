EXACT                   = 'yes'
EQUIVALENT              = 'yes_m'  # e.g. yes_modified
NO_MATCH                = 'no'


NULL_SET                = frozenset(['NULL','null','', '-', 'p.?', None])
SYN_ALIAS_SET           = frozenset(['p.(=)', 'p.='])
COSMIC_NULL             = 'p.?'
NON_CODING              = '-'
NULL                    = '-'
HGVS_SYN                = 'p.(=)'
DELINS                  = 'delins'
DEL                     = 'del'
INS                     = 'ins'
DUP                     = 'dup'
SUB                     = 'sub'
DELINS_SPLIT            = 'del_ins'  # for hgvs examples: c.34delTTinsAA
UNKNOWN                 = '?'
INVERSION               = 'inv'
SYNONYMOUS              = 'synonymous'
NONSENSE                = 'nonsense'
MISSENSE                = 'missense'
INFRAME                 = 'inframe'
FRAMESHIFT              = 'fs'
INFRAME_INSERTION       = 'inframe_insertion'
UPSTREAM                = 'upstream'
START_LOST              = 'start_lost'
STOP_LOST               = 'stop_lost'
DINUCLEOTIDE            = 'dinucleotide'
CANNOT_ASSESS           = '?'
EXTENSION               = 'ext'
INTRONIC                = 'intronic'
MULTI_SUB               = 'multi_sub'


DEL_SET = [DEL]
INS_SET = [DUP,INS]
SUB_SET = [NONSENSE, SYNONYMOUS, START_LOST, STOP_LOST]

CODING_START = 'c.'
PROTEIN_START = 'p.'
GENOMIC_START = 'g.'
CODING_START_LIST = ['r.','m.','c.','n.']

AMINO_ACID_SINGLETS = [
	'A', 'R', 'N', 'D', 'B',
	'C',
	'E', 'Q', 'Z', 'G',
	'H',
	'I',
	'L', 'K',
	'M',
	'F', 'P',
	'S',
	'T', 'W', 'Y',
	'V',
	'U',
	'X'
]

AMINO_ACID_TRIPLETS = [
	'Ala', 'Arg', 'Asn', 'Asp', 'Asx',
	'Cys',
	'Glu', 'Gln', 'Glx', 'Gly',
	'His',
	'Ile',
	'Leu', 'Lys',
	'Met',
	'Phe', 'Pro',
	'Ser',
	'Thr', 'Trp', 'Tyr',
	'Val',
	'Sec',
	'Ter'
]


INFRAME_SET             = frozenset(['inframe_variant', 'inframe_deletion', 'inframe_insertion'])
MISSENSE_SET            = frozenset(['non_synonymous_codon', 'missense_variant'])
UPSTREAM_SET            = frozenset(['upstream_variant', 'upstream_gene_variant', '2KB_upstream_variant', '5_prime_UTR_variant'])
SYN_SET                 = frozenset(['synonymous_codon', 'synonymous_variant'])
INFRAME_DELETION_SET    = frozenset(['disruptive_inframe_deletion', 'inframe_deletion'])
INFRAME_INSERTION_SET   = frozenset(['disruptive_inframe_insertion', 'inframe_insertion'])
START_SET               = frozenset(['initiator_codon_variant', 'start_lost'])
EFFECT_EQ_LIST          = [INFRAME_DELETION_SET, INFRAME_SET, UPSTREAM_SET, SYN_SET, MISSENSE_SET, START_SET]

