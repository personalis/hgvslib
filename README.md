compareHGVS
-----------------------------
compareHGVS is a package for parsing and comparing the equivalency of two 
variant strings described according to recommendations by the Human Genome 
Variation Society (HGVS).


Features
-----------------------------
This package provides the following functionalities:
   - Assess whether two HGVS strings describe the same variant
   - Normalize an HGVS string to a reduced form
   - Parse attributes of a variant object (e.g. location, transcript, bases,
     etc)

Current limitations
-----------------------------
Due to the complexity of HGVS syntax, this package currently does not support
comparisons of HGVS descriptions that require querying the accession sequence. 
For example, although ‘NM_006015.4:c.492_494dupCGC’ and 
'NM_006015.4:c.492_493insCGC’ describe the same variant, this package currently 
does not verify that CGC are duplicated bases according to the reference sequence. 

The package was developed for the purpose of comparing to a reference ‘truth’ 
test set, from which deviations in positional justification are assumed to be 
non-HGVS compliant.

If there are HGVS expressions that are not currently supported by this tool,
please make a comment on the github or send an email to
jennifer.yen@personalis.com.


Future plans
-------------------
-Inte

Unit Tests
-------------------

Over 115 unit tests have been created to test this code. 

To run:

$ python runtests.py -v


Example usage
-------------------

$ from hgvslib.class_functions import compare_hgvs

# the following check returns “yes”, indicating that this is an exact match.

$ hgvs1 = 'NM_000352.3:c.123delA'
$ hgvs2 = 'NM_000352.3:c.123delA'
$ compare_hgvs(hgvs1, hgvs2)
'yes'

# the following checks return “yes_m”, indicating that these are equivalent 
# but not exact matches ('yes modified').

$ hgvs1 = 'NP_005647.3:p.Q29fs'
$ hgvs2 = 'NP_005647.3:p.Q29Efs*10’
$ compare_hgvs(hgvs1, hgvs2)
'yes_m'

$ hgvs1 = 'p.Ser78Ser='
$ hgvs2 = 'p.='
$ compare_hgvs(hgvs1, hgvs2)
'yes_m'

$ hgvs1 = 'p.Glu78fs'
$ hgvs2 = 'p.Glu78GlyfsTer7'
$ compare_hgvs(hgvs1, hgvs2)
'yes_m'

$ hgvs1 = 'NM_004360.3:c.48+6_48+7delinsTT'
$ hgvs2 = 'NM_004360.3:c.48+6_48+7delCCinsTT'
$ compare_hgvs(hgvs1, hgvs2)
'yes_m'

# the following check returns “no”, indicating that these are different
# variants.

$ hgvs1 = 'NM_000352.3:c.123delA'
$ hgvs2 = 'NM_000352.3:c.125delA'
$ compare_hgvs(hgvs1, hgvs2)
'no'

$ hgvs1 = 'NP_005647.3:p.Q29fs'
$ hgvs2 = 'NP_005647.3:p.Q29Terfs'
$ compare_hgvs(hgvs1, hgvs2)
'no'

# the library can also perform basic parsing steps from an HGVS name.

$ from hgvslib.cHGVS import cHGVS

$ hgvs = cHGVS('NM_000352.3:c.215-10A>G')
$ hgvs.transcript
'NM_000352.3'
$ hgvs.name
'c.215-10A>G'
$ hgvs.type
'substitution'


Example script
-------------------

The example script reads in a text file with a reference HGVS string and 
columns of HGVS output from three different annotation callers. 
It compares the different HGVS strings in each row with the reference annotation
and outputs the result of the comparison: an exact match 'yes', equivalent 
'yes_m' or non-match 'no'.
 
$ python example/example_check_variants.py --infile example/test_file.txt



