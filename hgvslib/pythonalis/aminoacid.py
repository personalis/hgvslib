import os
import sys

import re

from hgvslib.pythonalis.mixins import Registry

AMINO_ACIDS = [
   'Alanine', 'Arginine', 'Asparagine', 'Aspartic_acid', 'Asparagine_Aspartic_Acid',
   'Cysteine',
   'Glutamic_Acid', 'Glutamine', 'Glutamine_Glutamic_Acid', 'Glycine',
   'Histidine',
   'Isoleucine',
   'Leucine', 'Lysine',
   'Methionine',
   'Phenylalanine', 'Proline',
   'Serine',
   'Threonine', 'Tryptophan', 'Tyrosine',
   'Valine',
   'Selenocysteine',
   'Termination'
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

class AminoAcid(Registry):
   Registry._add_registry('singlet_to_amino', {
      singlet: amino
      for amino, singlet
      in zip(AMINO_ACIDS, AMINO_ACID_SINGLETS)
   })

   Registry._add_registry('amino_to_singlet', {
      amino: singlet
      for amino, singlet
      in zip(AMINO_ACIDS, AMINO_ACID_SINGLETS)
   })

   Registry._add_registry('triplet_to_amino', {
      triplet: amino
      for amino, triplet
      in zip(AMINO_ACIDS, AMINO_ACID_TRIPLETS)
   })

   Registry._add_registry('amino_to_triplet', {
      amino: triplet
      for amino, triplet
      in zip(AMINO_ACIDS, AMINO_ACID_TRIPLETS)
   })

   @classmethod
   def singlet_from_amino(cls, amino_acid):
      return cls.registry_lookup('amino_to_singlet', amino_acid)

   @classmethod
   def amino_from_singlet(cls, amino_singlet):
      return cls.registry_lookup('singlet_to_amino', amino_singlet)

   @classmethod
   def triplet_from_amino(cls, amino_acid):
      return cls.registry_lookup('amino_to_triplet', amino_acid)

   @classmethod
   def amino_from_triplet(cls, amino_triplet):
      return cls.registry_lookup('triplet_to_amino', amino_triplet)

   @classmethod
   def triplet_from_singlet(cls, amino_singlet):
      return cls.triplet_from_amino(cls.amino_from_singlet(amino_singlet))

   @classmethod
   def singlet_from_triplet(cls, amino_triplet):
      # make sure the first letter is capitalized and the rest are lowercase.
      amino_triplet = amino_triplet[0].upper() + amino_triplet[1:].lower()
      return cls.singlet_from_amino(cls.amino_from_triplet(amino_triplet))