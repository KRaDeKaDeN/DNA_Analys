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
		"""Translate RNA to protein sequence"""
		genetic_code = {
			"AUG": "M",
			"UUU": "F",
			"UUC": "F",
			"UUA": "L",
			"UUG": "L",
			"UCU": "S",
			"UCC": "S",
			"UCA": "S",
			"UCG": "S",
			"UAU": "Y",
			"UAC": "Y",
			"UGU": "C",
			"UGC": "C",
			"UGG": "W",
			"CUU": "L",
			"CUC": "L",
			"CUA": "L",
			"CUG": "L",
			"CCU": "P",
			"CCC": "P",
			"CCA": "P",
			"CCG": "P",
			"CAU": "H",
			"CAC": "H",
			"CAA": "Q",
			"CAG": "Q",
			"CGU": "R",
			"CGC": "R",
			"CGA": "R",
			"CGG": "R",
			"AUU": "I",
			"AUC": "I",
			"AUA": "I",
			"GUU": "V",
			"GUC": "V",
			"GUA": "V",
			"GUG": "V",
			"ACU": "T",
			"ACC": "T",
			"ACA": "T",
			"ACG": "T",
			"AAU": "N",
			"AAC": "N",
			"AAA": "K",
			"AAG": "K",
			"GAU": "D",
			"GAC": "D",
			"GAA": "E",
			"GAG": "E",
			"GGU": "G",
			"GGC": "G",
			"GGA": "G",
			"GGG": "G",
			"UAA": "*",
			"UAG": "*",
			"UGA": "*",
		}
		protein = []
		for i in range(0, len(self) - 2, 3):
			codon = self.sequence[i : i + 3]
			amino_acid = genetic_code.get(codon, "?")
			if amino_acid == "*":
				break
			protein.append(amino_acid)
		return "".join(protein)
