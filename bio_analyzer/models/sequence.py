from abc import ABC, abstractmethod
from typing import Dict, Set


class SequenceValidationError(Exception):
	"""Base exception for sequence validation errors"""

	pass


class __Sequence(ABC):
	"""Abstract base class for biological sequences"""

	VALID_BASES: Set[str] = set()
	MIN_LENGTH: int = 3

	def __init__(self, sequence: str):
		self._sequence = sequence.upper()
		self.validate()

	def __len__(self) -> int:
		return len(self._sequence)

	@property
	def sequence(self) -> str:
		return self._sequence

	@sequence.setter
	def sequence(self, value: str):
		self._sequence = value.upper()
		self.validate()

	def validate(self):
		"""Validate sequence structure and content"""
		self._check_empty()
		self._check_length()
		self._check_invalid_chars()

	def _check_empty(self):
		if not self._sequence:
			raise SequenceValidationError("Sequence cannot be empty")

	def _check_length(self):
		if len(self) < self.MIN_LENGTH:
			raise SequenceValidationError(f"Minimum sequence length: {self.MIN_LENGTH}")

	def _check_invalid_chars(self):
		invalid = set(self._sequence) - self.VALID_BASES
		if invalid:
			raise SequenceValidationError(f"Invalid characters: {', '.join(invalid)}")

	@abstractmethod
	def count_bases(self) -> Dict[str, int]:
		"""Count individual bases in the sequence"""
		pass

	@abstractmethod
	def gc_content(self) -> float:
		"""Calculate GC content percentage"""
		pass


class DNA(__Sequence):
	"""DNA sequence implementation"""

	VALID_BASES = {"A", "T", "G", "C"}

	def count_bases(self) -> Dict[str, int]:
		return {base: self.sequence.count(base) for base in self.VALID_BASES}

	def gc_content(self) -> float:
		counts = self.count_bases()
		return ((counts["G"] + counts["C"]) / len(self)) * 100

	def transcribe(self) -> "RNA":
		"""Transcribe DNA to RNA"""
		return RNA(self.sequence.replace("T", "U"))

	def find_codons(self, codon_type: str) -> list[int]:
		"""Find positions of specific codon types"""
		codon_map = {"start": ["ATG"], "stop": ["TAA", "TAG", "TGA"]}
		targets = codon_map.get(codon_type.lower(), [])
		return [
			i for i in range(0, len(self) - 2, 3) if self.sequence[i : i + 3] in targets
		]

	def reverse_complement(self) -> str:
		"""Generate reverse complementary sequence"""
		complement = {"A": "T", "T": "A", "C": "G", "G": "C"}
		return "".join(complement.get(base, "") for base in reversed(self.sequence))

	def find_palindromic_sequences(
		self, min_len: int = 4, max_len: int = 12
	) -> list[dict]:
		results = []
		rev_comp = self.reverse_complement()

		# First we find all possible palindromes
		for length in range(max_len, min_len - 1, -1):  # Come from greater to smaller
			if length > len(self.sequence):
				continue

			for i in range(len(self.sequence) - length + 1):
				substr = self.sequence[i : i + length]
				if substr == rev_comp[-(i + length) : -i if i != 0 else None]:
					results.append(
						{
							"start": i,
							"end": i + length,
							"length": length,
							"sequence": substr,
						}
					)

		# Filtering of nested palindromes
		filtered = []
		results.sort(
			key=lambda x: (-x["length"], x["start"])
		)  # Sorting in length and position

		occupied = []
		for item in results:
			# We check if the current palindrome is intersected with the already selected
			overlap = any(
				item["start"] < existing["end"] and item["end"] > existing["start"]
				for existing in occupied
			)

			if not overlap:
				filtered.append(
					{
						"position": item["start"] + 1,
						"length": item["length"],
						"sequence": item["sequence"],
					}
				)
				occupied.append(item)

		return sorted(filtered, key=lambda x: x["position"])


class RNA(__Sequence):
	"""RNA sequence implementation"""

	VALID_BASES = {"A", "U", "G", "C"}

	def count_bases(self) -> Dict[str, int]:
		return {base: self.sequence.count(base) for base in self.VALID_BASES}

	def gc_content(self) -> float:
		counts = self.count_bases()
		return ((counts["G"] + counts["C"]) / len(self)) * 100

	def reverse_transcribe(self) -> DNA:
		"""Reverse transcribe RNA to DNA"""
		return DNA(self.sequence.replace("U", "T"))

	def translate(self) -> str:
		"""Translate RNA to protein sequence."""
		genetic_code = {
			"AUG": "M",  # Start codon
			"UUU": "F",
			"UUC": "F",  # Phenylalanine
			"UUA": "L",
			"UUG": "L",
			"CUU": "L",
			"CUC": "L",
			"CUA": "L",
			"CUG": "L",  # Leucine
			"UCU": "S",
			"UCC": "S",
			"UCA": "S",
			"UCG": "S",
			"AGU": "S",
			"AGC": "S",  # Serine
			"UAU": "Y",
			"UAC": "Y",  # Tyrosine
			"UGU": "C",
			"UGC": "C",  # Cysteine
			"UGG": "W",  # Tryptophan
			"CCU": "P",
			"CCC": "P",
			"CCA": "P",
			"CCG": "P",  # Proline
			"CAU": "H",
			"CAC": "H",  # Histidine
			"CAA": "Q",
			"CAG": "Q",  # Glutamine
			"CGU": "R",
			"CGC": "R",
			"CGA": "R",
			"CGG": "R",
			"AGA": "R",
			"AGG": "R",  # Arginine
			"AUU": "I",
			"AUC": "I",
			"AUA": "I",  # Isoleucine
			"GUU": "V",
			"GUC": "V",
			"GUA": "V",
			"GUG": "V",  # Valine
			"ACU": "T",
			"ACC": "T",
			"ACA": "T",
			"ACG": "T",  # Threonine
			"AAU": "N",
			"AAC": "N",  # Asparagine
			"AAA": "K",
			"AAG": "K",  # Lysine
			"GAU": "D",
			"GAC": "D",  # Aspartic acid
			"GAA": "E",
			"GAG": "E",  # Glutamic acid
			"GGU": "G",
			"GGC": "G",
			"GGA": "G",
			"GGG": "G",  # Glycine
			"GCU": "A",
			"GCC": "A",
			"GCA": "A",
			"GCG": "A",  # Alanine
			"UAA": "*",
			"UAG": "*",
			"UGA": "*",  # Stop codons
		}

		protein = []
		for i in range(0, len(self.sequence) - 2, 3):
			codon = self.sequence[i : i + 3]
			if len(codon) != 3:
				continue  # Skip incomplete codons
			amino_acid = genetic_code.get(codon, "?")
			protein.append(amino_acid)
		return "".join(protein)
