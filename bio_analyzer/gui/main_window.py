import tkinter as tk
import seaborn as sns
import pandas as pd
from tkinter import ttk, filedialog, messagebox
import matplotlib.pyplot as plt
from Bio import SeqIO
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
from ..services.analyzer import SequenceAnalyzer, AnalysisResult
from ..models.sequence import DNA, SequenceValidationError


class BioApp:
	"""Main application GUI"""

	def __init__(self, root):
		self.root = root
		self.root.title("Bio Sequence Analyzer")
		self.root.geometry("1200x800")
		self._current_result: AnalysisResult = None
		self.palindrome_min_len = tk.IntVar(value=4)
		self.palindrome_max_len = tk.IntVar(value=12)

		self._setup_ui()

	def _setup_ui(self):
		"""Initialize UI components"""
		self.notebook = ttk.Notebook(self.root)

		# Input Tab
		self.input_frame = ttk.Frame(self.notebook)
		self._create_input_tab()

		# Results Tab
		self.results_frame = ttk.Frame(self.notebook)
		self._create_results_tab()

		# Visualization Tab
		self.viz_frame = ttk.Frame(self.notebook)
		self._create_visualization_tab()

		self.palindrome_tab = ttk.Frame(self.notebook)
		self._create_palindrome_components()

		self.notebook.add(self.input_frame, text="Input")
		self.notebook.add(self.results_frame, text="Results")
		self.notebook.add(self.viz_frame, text="Visualization")
		self.notebook.add(self.palindrome_tab, text="Palindromes")
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

	def _create_palindrome_components(self):
		"""Visualization tab widgets"""
		main_frame = ttk.Frame(self.palindrome_tab)
		main_frame.pack(pady=20, padx=20, fill=tk.BOTH, expand=True)

		# Fields input parameters
		params_frame = ttk.Frame(main_frame)
		params_frame.pack(fill=tk.X, pady=10)

		ttk.Label(params_frame, text="Min length:").pack(side=tk.LEFT)
		min_entry = ttk.Entry(
			params_frame, textvariable=self.palindrome_min_len, width=5
		)
		min_entry.pack(side=tk.LEFT, padx=5)

		ttk.Label(params_frame, text="Max length:").pack(side=tk.LEFT)
		max_entry = ttk.Entry(
			params_frame, textvariable=self.palindrome_max_len, width=5
		)
		max_entry.pack(side=tk.LEFT, padx=5)

		# The area of ​​the results
		self.palindrome_result = ttk.Frame(main_frame)
		self.palindrome_result.pack(fill=tk.BOTH, expand=True)

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
		try:
			# Getting parameters from the user
			min_len = self.palindrome_min_len.get()
			max_len = self.palindrome_max_len.get()

			# Validation of parameters
			if min_len < 4:
				raise ValueError("Minimum length must be at least 4")
			if max_len < min_len:
				raise ValueError("Max length must be greater than min length")

			dna = DNA(self.sequence_entry.get(1.0, tk.END).strip())
			analyzer = SequenceAnalyzer(dna)
			analyzer.palindrome_params = (min_len, max_len)
			self._current_result = analyzer.analyze()
			self._show_results()
			self._update_visualizations()
			self._update_palindrome_display()

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
		self._create_codon_heatmap()
		# self._create_codon_chart()

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

	def _create_codon_heatmap(self):
		rna_seq = self._current_result.data.get("rna", "")
		if not rna_seq or len(rna_seq) < 3:
			return

		# Dictionary of amino acids and corresponding codons (RNA)
		amino_acid_groups = {
			"Ala (A)": ["GCU", "GCC", "GCA", "GCG"],
			"Arg (R)": ["CGU", "CGC", "CGA", "CGG", "AGA", "AGG"],
			"Asn (N)": ["AAU", "AAC"],
			"Asp (D)": ["GAU", "GAC"],
			"Cys (C)": ["UGU", "UGC"],
			"Gln (Q)": ["CAA", "CAG"],
			"Glu (E)": ["GAA", "GAG"],
			"Gly (G)": ["GGU", "GGC", "GGA", "GGG"],
			"His (H)": ["CAU", "CAC"],
			"Ile (I)": ["AUU", "AUC", "AUA"],
			"Leu (L)": ["UUA", "UUG", "CUU", "CUC", "CUA", "CUG"],
			"Lys (K)": ["AAA", "AAG"],
			"Met (M)": ["AUG"],
			"Phe (F)": ["UUU", "UUC"],
			"Pro (P)": ["CCU", "CCC", "CCA", "CCG"],
			"Ser (S)": ["UCU", "UCC", "UCA", "UCG", "AGU", "AGC"],
			"Thr (T)": ["ACU", "ACC", "ACA", "ACG"],
			"Trp (W)": ["UGG"],
			"Tyr (Y)": ["UAU", "UAC"],
			"Val (V)": ["GUU", "GUC", "GUA", "GUG"],
			"Stop (*)": ["UAA", "UAG", "UGA"],
		}

		# Calculation of the frequency of codons
		codon_counts = {}
		total_codons = len(rna_seq) // 3
		for i in range(total_codons):
			codon = rna_seq[i * 3 : (i + 1) * 3].upper()
			codon_counts[codon] = codon_counts.get(codon, 0) + 1

		# Creating Dataframe
		data = []
		for aa, codons in amino_acid_groups.items():
			for codon in codons:
				data.append(
					{
						"Amino Acid": aa,
						"Codon": codon,
						"Frequency": codon_counts.get(codon, 0),
					}
				)

		df = pd.DataFrame(data)

		# Creating a figure
		fig = plt.Figure(figsize=(13, 8), dpi=100)
		ax = fig.add_subplot(111)

		# Creating a pivot table
		pivot_df = df.pivot_table(
			index="Amino Acid",
			columns="Codon",
			values="Frequency",
			fill_value=0,
			aggfunc="sum",
		)

		# Sorting amino acids by frequency of use
		sorted_aa = pivot_df.sum(axis=1).sort_values(ascending=False).index
		pivot_df = pivot_df.loc[sorted_aa]

		# Sorting codons inside amino acids
		for aa in pivot_df.index:
			sorted_codons = pivot_df.loc[aa].sort_values(ascending=False).index
			pivot_df.loc[aa] = pivot_df.loc[aa][sorted_codons]

		# Visualization
		sns.heatmap(
			pivot_df,
			cmap="YlGnBu",
			annot=True,
			fmt="d",
			linewidths=0.4,
			ax=ax,
			cbar_kws={"label": "Frequency of Usage"},
		)

		# Office Settings
		ax.set_title("Codon Usage Heatmap (Grouped by Amino Acid)", pad=20, fontsize=14)
		ax.set_xlabel("Codon", fontsize=12)
		ax.set_ylabel("Amino Acid", fontsize=12)
		plt.xticks(rotation=45, ha="right", fontsize=8)
		plt.yticks(rotation=0, fontsize=8)

		# Integration in Gui
		canvas = FigureCanvasTkAgg(fig, master=self.viz_container)
		canvas.draw()
		canvas.get_tk_widget().pack(side=tk.TOP, fill=tk.BOTH, expand=True)

		# Processing of size changes
		def on_resize(event):
			fig.tight_layout()
			canvas.draw()

		canvas.get_tk_widget().bind("<Configure>", on_resize)

	def _update_palindrome_display(self):
		# Cleaning previous widgets
		for widget in self.palindrome_result.winfo_children():
			widget.destroy()

		# Verification of the availability of data
		if not self._current_result or not self._current_result._data["palindromes"]:
			no_data_label = ttk.Label(
				self.palindrome_result,
				text="No palindromic sequences found for current parameters",
				font=("Arial", 10, "italic"),
				foreground="gray50",
			)
			no_data_label.pack(pady=40, expand=True)
			return

		# Creating a container for a table and scrolling
		container = ttk.Frame(self.palindrome_result)
		container.pack(fill=tk.BOTH, expand=True, padx=10, pady=10)

		# Creating TreeView with vertical scrolling
		tree_frame = ttk.Frame(container)
		tree_frame.pack(fill=tk.BOTH, expand=True)

		columns = ("position", "length", "sequence")
		tree = ttk.Treeview(
			tree_frame, columns=columns, show="headings", selectmode="browse", height=12
		)

		# Configuration of the columns
		tree.heading("position", text="Position (1-based)", anchor=tk.CENTER)
		tree.heading("length", text="Length", anchor=tk.CENTER)
		tree.heading("sequence", text="Sequence", anchor=tk.W)

		tree.column("position", width=140, anchor=tk.CENTER, stretch=False)
		tree.column("length", width=80, anchor=tk.CENTER, stretch=False)
		tree.column("sequence", anchor=tk.W, stretch=True)

		# Setting styles for even/odd lines
		tree.tag_configure("even", background="#f8f8f8")
		tree.tag_configure("odd", background="white")

		# Adding data
		for idx, palindrome in enumerate(self._current_result._data["palindromes"]):
			tags = ("even",) if idx % 2 == 0 else ("odd",)
			tree.insert(
				"",
				tk.END,
				values=(
					palindrome["position"],
					palindrome["length"],
					palindrome["sequence"],
				),
				tags=tags,
			)

		# Scrolling Settings
		scrollbar = ttk.Scrollbar(tree_frame, orient=tk.VERTICAL, command=tree.yview)
		tree.configure(yscrollcommand=scrollbar.set)

		# Packaging of elements
		tree.pack(side=tk.LEFT, fill=tk.BOTH, expand=True)
		scrollbar.pack(side=tk.RIGHT, fill=tk.Y)

		# Adding an explanatory text
		info_label = ttk.Label(
			container,
			text="* Only longest non-overlapping palindromes are shown",
			foreground="gray40",
			font=("Arial", 8),
		)
		info_label.pack(side=tk.BOTTOM, pady=(5, 0))

		# Blocking Editing
		tree.bind("<Button-1>", lambda e: "break")
