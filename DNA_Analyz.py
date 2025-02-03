import tkinter as tk
from tkinter import ttk, filedialog, messagebox
from Bio import SeqIO
import matplotlib.pyplot as plt
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg


class DNAAnalyzer:

    def __init__(self, dna_sourse):
        self.dna_sourse = dna_sourse.upper()

    def validate_sequence(self):
        """Check that sequence contains only ATGC characters"""
        if not self.dna_sourse:
            raise ValueError("The sequence is empty.")
        invalid_chars = set(self.dna_sourse) - set("ATGC")
        if invalid_chars:
            raise ValueError(f"Invalid characters found: {invalid_chars}")
        return True

    def count_nucleotides(self):
        '''Calculates and returns the number of nucleotides in the sequence'''
        return {
            "A": self.dna_sourse.count('A'),
            "T": self.dna_sourse.count('T'),
            "G": self.dna_sourse.count('G'),
            "C": self.dna_sourse.count('C')
        }

    def gc_content(self):
        '''Calculates the percentage of the content of nucletide pair GC'''
        counts = self.count_nucleotides()
        gc_count = counts["G"] + counts["C"]
        return (gc_count / len(self.dna_sourse) * 100)

    def transcribe(self):
        '''Transcribes the introduced DNA in RNA replacing Timin with Uracyl'''
        return self.dna_sourse.replace('T', 'U')

    def find_codons(self, codons_type='start'):
        '''Finds start or stop codon'''
        codons = []
        if codons_type == 'start':
            target_codons = ['ATG']
        elif codons_type == 'stop':
            target_codons = ['TAA', 'TAG', 'TGA']
        else:
            return []
        for i in range(0, len(self.dna_sourse) - 2, 3):
            codon = self.dna_sourse[i:i + 3]
            if codon in target_codons:
                codons.append(i)
        return codons

    @staticmethod
    def translate_rna_to_protein(rna_sequence):
        '''Translate RNA into protein, genetic_code contains the accordance of triplets to their aminoacids'''
        genetic_code = {
            'AUG': 'M',
            'UUU': 'F',
            'UUC': 'F',
            'UUA': 'L',
            'UUG': 'L',
            'UCU': 'S',
            'UCC': 'S',
            'UCA': 'S',
            'UCG': 'S',
            'UAU': 'Y',
            'UAC': 'Y',
            'UGU': 'C',
            'UGC': 'C',
            'UGG': 'W',
            'CUU': 'L',
            'CUC': 'L',
            'CUA': 'L',
            'CUG': 'L',
            'CCU': 'P',
            'CCC': 'P',
            'CCA': 'P',
            'CCG': 'P',
            'CAU': 'H',
            'CAC': 'H',
            'CAA': 'Q',
            'CAG': 'Q',
            'CGU': 'R',
            'CGC': 'R',
            'CGA': 'R',
            'CGG': 'R',
            'AUU': 'I',
            'AUC': 'I',
            'AUA': 'I',
            'GUU': 'V',
            'GUC': 'V',
            'GUA': 'V',
            'GUG': 'V',
            'ACU': 'T',
            'ACC': 'T',
            'ACA': 'T',
            'ACG': 'T',
            'AAU': 'N',
            'AAC': 'N',
            'AAA': 'K',
            'AAG': 'K',
            'GAU': 'D',
            'GAC': 'D',
            'GAA': 'E',
            'GAG': 'E',
            'GGU': 'G',
            'GGC': 'G',
            'GGA': 'G',
            'GGG': 'G',
            'UAA': '*',
            'UAG': '*',
            'UGA': '*'
        }

        protein = []
        for i in range(0, len(rna_sequence) - 2, 3):
            codon = rna_sequence[i:i + 3]
            if codon in genetic_code:
                amino_acid = genetic_code[codon]
                if amino_acid == '*':
                    break
                protein.append(amino_acid)
        else:
            protein.append('?')
        return ''.join(protein)


class DNAApp:

    def __init__(self, root):
        self.root = root
        self.root.title("DNA Analysis Tool")
        self.dna_source = ""

        # Create tab control
        self.notebook = ttk.Notebook(root)
        self.notebook.pack(fill='both', expand=True)

        # Create tabs
        self.create_input_tab()
        self.create_analysis_tab()
        self.create_plots_tab()

    def create_input_tab(self):
        """Create data input tab"""
        self.input_tab = ttk.Frame(self.notebook)
        self.notebook.add(self.input_tab, text="Data Input")

        # Sequence input field
        tk.Label(self.input_tab, text="Enter DNA sequence:").pack(pady=5)
        self.sequence_text = tk.Text(self.input_tab, height=10, width=80)
        self.sequence_text.pack(pady=5)

        # File load button
        ttk.Button(self.input_tab, text="Load FASTA",
                   command=self.load_fasta).pack(pady=5)

    def create_analysis_tab(self):
        """Create analysis results tab"""
        self.analysis_tab = ttk.Frame(self.notebook)
        self.notebook.add(self.analysis_tab, text="Analysis")

        # Results display area
        self.results_text = tk.Text(self.analysis_tab, height=15, width=80)
        self.results_text.pack(pady=10, padx=10)

        # Analysis button
        ttk.Button(self.analysis_tab,
                   text="Run Analysis",
                   command=self.run_analysis).pack(pady=5)

    def create_plots_tab(self):
        """Create visualization tab"""
        self.plots_tab = ttk.Frame(self.notebook)
        self.notebook.add(self.plots_tab, text="Visualizations")

        # Plot frames
        self.plot_frame1 = ttk.Frame(self.plots_tab)
        self.plot_frame1.pack(side='top', fill='both', expand=True)

        self.plot_frame2 = ttk.Frame(self.plots_tab)
        self.plot_frame2.pack(side='top', fill='both', expand=True)

    def load_fasta(self):
        """Load FASTA file"""
        filepath = filedialog.askopenfilename(filetypes=[("FASTA files",
                                                          "*.fasta *.fa")])
        if not filepath: 
            return

        try:
            with open(filepath) as f:
                record = next(SeqIO.parse(f, "fasta"))
                self.dna_source = str(record.seq)
                self.sequence_text.delete(1.0, tk.END)
                self.sequence_text.insert(tk.END, self.dna_source)
        except Exception as e:
            messagebox.showerror("Error", f"Error loading file:\n{str(e)}")

    def run_analysis(self):
        """Run sequence analysis"""
        self.dna_source = self.sequence_text.get(1.0, tk.END).strip().upper()

        try:
            analyzer = DNAAnalyzer(self.dna_source)
            analyzer.validate_sequence()

            # Perform calculations
            counts = analyzer.count_nucleotides()
            gc = analyzer.gc_content()
            rna = analyzer.transcribe()
            protein = analyzer.translate_rna_to_protein(rna)

            # Generate report
            report = (f"Sequence length: {len(self.dna_source)}\n"
                      f"Nucleotide composition:\n"
                      f"A: {counts['A']} | T: {counts['T']}\n"
                      f"G: {counts['G']} | C: {counts['C']}\n"
                      f"GC content: {gc:.2f}%\n\n"
                      f"Start codons: {analyzer.find_codons('start')}\n"
                      f"Stop codons: {analyzer.find_codons('stop')}\n\n"
                      f"RNA: {rna[:50]}...\n\n"
                      f"Protein: {protein}")

            self.results_text.delete(1.0, tk.END)
            self.results_text.insert(tk.END, report)
            self.update_plots(analyzer)

        except ValueError as e:
            messagebox.showerror("Error", str(e))

    def update_plots(self, analyzer):
        """Update visualization plots"""
        # Clear previous plots
        for frame in [self.plot_frame1, self.plot_frame2]:
            for widget in frame.winfo_children():
                widget.destroy()

        # Nucleotide distribution plot
        counts = analyzer.count_nucleotides()
        fig1 = plt.Figure(figsize=(4, 4))
        ax1 = fig1.add_subplot(111)
        ax1.pie(counts.values(), labels=counts.keys(), autopct='%1.1f%%')
        ax1.set_title("Nucleotide Distribution")
        canvas1 = FigureCanvasTkAgg(fig1, self.plot_frame1)
        canvas1.get_tk_widget().pack(side='left', fill='both', expand=True)

        # Codon usage frequency plot
        rna = analyzer.transcribe()
        codon_counts = {}
        for i in range(0, len(rna) - 2, 3):
            codon = rna[i:i + 3]
            codon_counts[codon] = codon_counts.get(codon, 0) + 1

        fig2 = plt.Figure(figsize=(8, 4))
        ax2 = fig2.add_subplot(111)
        ax2.bar(codon_counts.keys(), codon_counts.values())
        ax2.set_title("Codon Usage Frequency")
        ax2.tick_params(axis='x', rotation=45)
        canvas2 = FigureCanvasTkAgg(fig2, self.plot_frame2)
        canvas2.get_tk_widget().pack(fill='both', expand=True)


if __name__ == "__main__":
    root = tk.Tk()
    app = DNAApp(root)
    root.geometry("900x900")
    root.mainloop()
