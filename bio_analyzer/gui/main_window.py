import tkinter as tk
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
        self.current_result: AnalysisResult = None
        
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
        
        self.notebook.add(self.input_frame, text="Input")
        self.notebook.add(self.results_frame, text="Results")
        self.notebook.add(self.viz_frame, text="Visualization")
        self.notebook.pack(expand=True, fill='both')
    def _create_input_tab(self):
        """Creating the input tab"""
        frame = ttk.Frame(self.input_frame)
        frame.pack(pady=20, padx=20, fill='both', expand=True)

        ttk.Label(frame, text="DNA Sequence Input:").pack(anchor='w')
        self.sequence_entry = tk.Text(frame, height=10, width=100)
        self.sequence_entry.pack(pady=10)

        btn_frame = ttk.Frame(frame)
        btn_frame.pack(fill='x')

        ttk.Button(btn_frame, text="Load FASTA",
                   command=self._load_fasta).pack(side='left')
        ttk.Button(btn_frame, text="Analyze",
                   command=self._run_analysis).pack(side='right')

    def _create_results_tab(self):
        """creating results tab"""
        frame = ttk.Frame(self.results_frame)
        frame.pack(fill='both', expand=True, padx=20, pady=20)

        self.results_text = tk.Text(frame, wrap=tk.WORD)
        scrollbar = ttk.Scrollbar(frame, command=self.results_text.yview)
        self.results_text.configure(yscrollcommand=scrollbar.set)

        self.results_text.pack(side='left', fill='both', expand=True)
        scrollbar.pack(side='right', fill='y')

    def _create_visualization_tab(self):
        """Creating a visualization tab"""
        self.viz_container = ttk.Frame(self.viz_frame)
        self.viz_container.pack(fill='both', expand=True, padx=20, pady=20)

    def _load_fasta(self):
        """Upload FASTA file"""
        filepath = filedialog.askopenfilename(
            filetypes=[("FASTA files", "*.fasta *.fa"), ("All files", "*.*")])
        if not filepath:
            return

        try:
            with open(filepath) as f:
                record = next(SeqIO.parse(f, "fasta"))
                self.sequence_entry.delete(1.0, tk.END)
                self.sequence_entry.insert(tk.END, str(record.seq))
                messagebox.showinfo("Success",
                                    "FASTA file loaded successfully!")
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
            messagebox.showerror("Analysis Error",
                                 f"Unexpected error: {str(e)}")

    def _show_results(self):
        """Display of results"""
        if self._current_result:
            self.results_text.delete(1.0, tk.END)
            self.results_text.insert(tk.END,
                                     self._current_result.format_report())

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
        counts = self._current_result.data['base_counts']
        fig = plt.Figure(figsize=(5, 4))
        ax = fig.add_subplot(111)
        ax.pie(counts.values(),
               labels=counts.keys(),
               autopct='%1.1f%%',
               colors=['#ff9999', '#66b3ff', '#99ff99', '#ffcc99'])
        ax.set_title("Nucleotide Distribution")

        canvas = FigureCanvasTkAgg(fig, self.viz_container)
        canvas.get_tk_widget().pack(side='left', fill='both', expand=True)

    def _create_codon_chart(self):
        """Creating a histogram of codons"""
        rna_seq = self._current_result.data.get('rna')
        if not rna_seq:
            return

        codon_counts = {}
        for i in range(0, len(rna_seq) - 2, 3):
            codon = rna_seq[i:i + 3]
            codon_counts[codon] = codon_counts.get(codon, 0) + 1

        fig = plt.Figure(figsize=(10, 4))
        ax = fig.add_subplot(111)
        ax.bar(codon_counts.keys(), codon_counts.values(), color='skyblue')
        ax.set_title("Codon Usage Frequency")
        ax.tick_params(axis='x', rotation=45)

        canvas = FigureCanvasTkAgg(fig, self.viz_container)
        canvas.get_tk_widget().pack(side='left', fill='both', expand=True)