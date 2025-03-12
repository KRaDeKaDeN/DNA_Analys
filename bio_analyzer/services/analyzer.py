from typing import Dict
from ..models.sequence import DNA


class AnalysisResult:
	"""Container for analysis results"""

	def __init__(self, data: Dict):
		self._data = data

	@property
	def data(self) -> Dict:
		return self._data

	def format_report(self) -> str:
		"""Format results for display"""
		return (
			f"DNA Sequence Analysis Report\n{'=' * 40}\n"
			f"Sequence Length: {self._data['length']}\n\n"
			f"Nucleotide Composition:\n"
			f"A: {self._data['base_counts']['A']} | T: {self._data['base_counts']['T']}\n"
			f"G: {self._data['base_counts']['G']} | C: {self._data['base_counts']['C']}\n"
			f"GC Content: {self._data['gc_content']:.2f}%\n\n"
			f"Start Codons Positions: {self._data['start_codons']}\n"
			f"Stop Codons Positions: {self._data['stop_codons']}\n\n"
			f"RNA Sequence (first 50 bases):\n{self._data['rna'][:50] if self._data['rna'] else 'N/A'}...\n\n"
			f"Protein Sequence:\n{self._data['protein'] or 'N/A'}\n\n"
			f"Palindromic Sequences:\n {len(self._data['palindromes'])} found\n"
		)


class SequenceAnalyzer:
	"""Main analysis service"""

	def __init__(self, dna: DNA):
		self._dna = dna
		self._rna = dna.transcribe()
		self._protein = self._rna.translate()

	def set_palindrome_params(self, min_len: int, max_len: int):
		"""Checking and installing parameters"""
		if min_len < 4:
			raise ValueError("Minimum length must be ≥4")
		if max_len < min_len:
			raise ValueError("Max length must be ≥ min length")
		if max_len > 1000:
			raise ValueError("For sequences >1000bp use batch processing")
		self.palindrome_params = (min_len, max_len)

	def analyze(self) -> AnalysisResult:
		"""Perform comprehensive analysi swith current parameters"""
		min_len, max_len = self.palindrome_params
		return AnalysisResult(
			{
				"length": len(self._dna),
				"base_counts": self._dna.count_bases(),
				"gc_content": self._dna.gc_content(),
				"start_codons": self._dna.find_codons("start"),
				"stop_codons": self._dna.find_codons("stop"),
				"rna": self._rna.sequence,
				"protein": self._protein,
				"palindromes": self._dna.find_palindromic_sequences(min_len, max_len),
			}
		)
