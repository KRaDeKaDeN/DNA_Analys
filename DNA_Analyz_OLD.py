import tkinter as tk
from tkinter import ttk, filedialog, messagebox
from abc import ABC, abstractmethod
from Bio import SeqIO
import matplotlib.pyplot as plt
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
from typing import Dict, List, Optional


class SequenceValidationError(Exception):
	"""The basic exception for the errors of the validation of sequences"""

	pass


class Sequence(ABC):
	"""Abstract basic class for biological sequences"""

	VALID_BASES: set = set()
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
		"""sequence validation"""
		self._check_empty()
		self._check_length()
		self._check_invalid_chars()

	def _check_empty(self):
		if not self._sequence:
			raise SequenceValidationError("Sequence cannot be empty")

	def _check_length(self):
		if len(self) < self.MIN_LENGTH:
			raise SequenceValidationError(
				f"Sequence too short. Minimum length: {self.MIN_LENGTH}"
			)

	def _check_invalid_chars(self):
		invalid = set(self._sequence) - self.VALID_BASES
		if invalid:
			raise SequenceValidationError(f"Invalid characters: {', '.join(invalid)}")

	@abstractmethod
	def count_bases(self) -> Dict[str, int]:
		"""nucleotides count"""
		pass

	@abstractmethod
	def gc_content(self) -> float:
		"""GC-content"""
		pass


class DNA(Sequence):
	"""Class for working with DNA sequences"""

	VALID_BASES = {"A", "T", "G", "C"}

	def count_bases(self) -> Dict[str, int]:
		return {base: self.sequence.count(base) for base in self.VALID_BASES}

	def gc_content(self) -> float:
		counts = self.count_bases()
		return ((counts["G"] + counts["C"]) / len(self)) * 100

	def transcribe(self) -> "RNA":
		"""DNA transcription to RNA"""
		return RNA(self.sequence.replace("T", "U"))

	def find_codons(self, codon_type: str) -> List[int]:
		"""Search for codons of a given type"""
		codon_map = {"start": ["ATG"], "stop": ["TAA", "TAG", "TGA"]}
		targets = codon_map.get(codon_type.lower(), [])
		return [
			i for i in range(0, len(self) - 2, 3) if self.sequence[i : i + 3] in targets
		]


class RNA(Sequence):
	"""Class for working with RNA sequences"""

	VALID_BASES = {"A", "U", "G", "C"}

	def count_bases(self) -> Dict[str, int]:
		return {base: self.sequence.count(base) for base in self.VALID_BASES}

	def gc_content(self) -> float:
		counts = self.count_bases()
		return ((counts["G"] + counts["C"]) / len(self)) * 100

	def reverse_transcribe(self) -> DNA:
		"""Обратная транскрипция РНК в ДНК"""
		return DNA(self.sequence.replace("U", "T"))

	def translate(self) -> str:
		"""RNA translation to protein"""
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


# --- Модуль анализа данных ---
class AnalysisResult:
	"""Class for storing analysis results"""

	def __init__(self, data: Dict):
		self._data = data

	@property
	def data(self) -> Dict:
		return self._data

	def format_report(self) -> str:
		"""Formatting a report for output"""
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
			f"Protein Sequence:\n{self._data['protein'] or 'N/A'}"
		)


class SequenceAnalyzer:
	"""Class for performing sequences"""

	def __init__(self, dna: DNA):
		self._dna = dna
		self._rna = dna.transcribe()
		self._protein = self._rna.translate()

	def analyze(self) -> AnalysisResult:
		"""The main method of analysis"""
		return AnalysisResult(
			{
				"length": len(self._dna),
				"base_counts": self._dna.count_bases(),
				"gc_content": self._dna.gc_content(),
				"start_codons": self._dna.find_codons("start"),
				"stop_codons": self._dna.find_codons("stop"),
				"rna": self._rna.sequence,
				"protein": self._protein,
			}
		)


class BioApp:
	"""Graphic application interface"""

	def __init__(self, root):
		self.root = root
		self.root.title("Advanced Bio Analyzer")
		self.root.geometry("1200x800")

		self._setup_ui()
		self._current_result: Optional[AnalysisResult] = None

	def _setup_ui(self):
		"""interface initialisation"""
		self.notebook = ttk.Notebook(self.root)

		# Вкладка ввода
		self.input_frame = ttk.Frame(self.notebook)
		self._create_input_tab()

		# Вкладка результатов
		self.results_frame = ttk.Frame(self.notebook)
		self._create_results_tab()

		# Вкладка визуализации
		self.viz_frame = ttk.Frame(self.notebook)
		self._create_visualization_tab()

		self.notebook.add(self.input_frame, text="Input")
		self.notebook.add(self.results_frame, text="Results")
		self.notebook.add(self.viz_frame, text="Visualization")
		self.notebook.pack(expand=True, fill="both")

	def _create_input_tab(self):
		"""Creating the input tab"""
		frame = ttk.Frame(self.input_frame)
		frame.pack(pady=20, padx=20, fill="both", expand=True)

		ttk.Label(frame, text="DNA Sequence Input:").pack(anchor="w")
		self.sequence_entry = tk.Text(frame, height=10, width=100)
		self.sequence_entry.pack(pady=10)

		btn_frame = ttk.Frame(frame)
		btn_frame.pack(fill="x")

		ttk.Button(btn_frame, text="Load FASTA", command=self._load_fasta).pack(
			side="left"
		)
		ttk.Button(btn_frame, text="Analyze", command=self._run_analysis).pack(
			side="right"
		)

	def _create_results_tab(self):
		"""creating results tab"""
		frame = ttk.Frame(self.results_frame)
		frame.pack(fill="both", expand=True, padx=20, pady=20)

		self.results_text = tk.Text(frame, wrap=tk.WORD)
		scrollbar = ttk.Scrollbar(frame, command=self.results_text.yview)
		self.results_text.configure(yscrollcommand=scrollbar.set)

		self.results_text.pack(side="left", fill="both", expand=True)
		scrollbar.pack(side="right", fill="y")

	def _create_visualization_tab(self):
		"""Creating a visualization tab"""
		self.viz_container = ttk.Frame(self.viz_frame)
		self.viz_container.pack(fill="both", expand=True, padx=20, pady=20)

	def _load_fasta(self):
		"""Upload FASTA file"""
		filepath = filedialog.askopenfilename(
			filetypes=[("FASTA files", "*.fasta *.fa"), ("All files", "*.*")]
		)
		if not filepath:
			return

		try:
			with open(filepath) as f:
				record = next(SeqIO.parse(f, "fasta"))
				self.sequence_entry.delete(1.0, tk.END)
				self.sequence_entry.insert(tk.END, str(record.seq))
				messagebox.showinfo("Success", "FASTA file loaded successfully!")
		except Exception as e:
			messagebox.showerror("Error", f"Error loading file:\n{str(e)}")

	def _run_analysis(self):
		"""Analysis"""
		try:
			dna = DNA(self.sequence_entry.get(1.0, tk.END).strip())
			analyzer = SequenceAnalyzer(dna)
			self._current_result = analyzer.analyze()
			self._show_results()
			self._update_visualizations()
		except SequenceValidationError as e:
			messagebox.showerror("Validation Error", str(e))
		except Exception as e:
			messagebox.showerror("Analysis Error", f"Unexpected error: {str(e)}")

	def _show_results(self):
		"""Display of results"""
		if self._current_result:
			self.results_text.delete(1.0, tk.END)
			self.results_text.insert(tk.END, self._current_result.format_report())

	def _update_visualizations(self):
		"""Updating visualizations"""
		# Cleaning previous graphs
		for widget in self.viz_container.winfo_children():
			widget.destroy()

		if not self._current_result:
			return

		# Creation of new graphs
		self._create_pie_chart()
		self._create_codon_chart()

	def _create_pie_chart(self):
		"""Creating a circular diagram"""
		counts = self._current_result.data["base_counts"]
		fig = plt.Figure(figsize=(5, 4))
		ax = fig.add_subplot(111)
		ax.pie(
			counts.values(),
			labels=counts.keys(),
			autopct="%1.1f%%",
			colors=["#ff9999", "#66b3ff", "#99ff99", "#ffcc99"],
		)
		ax.set_title("Nucleotide Distribution")

		canvas = FigureCanvasTkAgg(fig, self.viz_container)
		canvas.get_tk_widget().pack(side="left", fill="both", expand=True)

	def _create_codon_chart(self):
		"""Creating a histogram of codons"""
		rna_seq = self._current_result.data.get("rna")
		if not rna_seq:
			return

		codon_counts = {}
		for i in range(0, len(rna_seq) - 2, 3):
			codon = rna_seq[i : i + 3]
			codon_counts[codon] = codon_counts.get(codon, 0) + 1

		fig = plt.Figure(figsize=(10, 4))
		ax = fig.add_subplot(111)
		ax.bar(codon_counts.keys(), codon_counts.values(), color="skyblue")
		ax.set_title("Codon Usage Frequency")
		ax.tick_params(axis="x", rotation=45)

		canvas = FigureCanvasTkAgg(fig, self.viz_container)
		canvas.get_tk_widget().pack(side="left", fill="both", expand=True)


if __name__ == "__main__":
	root = tk.Tk()
	app = BioApp(root)
	root.mainloop()
