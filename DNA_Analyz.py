import tkinter as tk
from tkinter import filedialog, messagebox
from Bio import SeqIO
import matplotlib.pyplot as plt
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg


class DNAAnalyzer:

    def __init__(self, dna_sourse):
        self.dna_sourse = dna_sourse.upper()

    def validate_sequence(self):
        return all(base in "ATGC" for base in self.dna_sourse)

    def count_nucleotides(self):
        return {
            "A": self.dna_sourse.count('A'),
            "T": self.dna_sourse.count('T'),
            "G": self.dna_sourse.count('G'),
            "C": self.dna_sourse.count('C')
        }

    def gc_content(self):
        counts = self.count_nucleotides()
        gc_count = counts["G"] + counts["C"]
        return (gc_count / len(self.dna_sourse))

    def transcribe(self):
        return self.dna_sourse.replace('T', 'U')

    def find_codons(self, codons_type='start'):
        codons = []
        if codons_type == 'start':
            target_codons = ['ATG']
        elif codons_type == 'stop':
            target_codons = ['TAA', 'TAG', 'TGA']
        else:
            return []
        for i in range(0, len(self.dna_sourse) - 2, 3):
            codon = self.dna_sourse[i:i + 3]
            if codon in target_codons: codons.append(i)
        return codons

    @staticmethod
    def translate_rna_to_protein(rna_sequence):
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
        self.root.title("Analyzis DNA and RNA")
        self.root.geometry("800x600")
        self.dna_sourse = ""
        self.result_text = tk.StringVar()
        self.create_widgets()

    def create_widgets(self):

        label_instruction = tk.Label(
            self.root, text="Enter DNA sequence or upload FASTA file:")
        label_instruction.pack(pady=10)
        self.entry_sequence = tk.Entry(self.root, width=80)
        self.entry_sequence.pack(pady=5)
        btn_analyse = tk.Button(self.root,
                                text="Analyze",
                                command=self.analyze_sequence)
        btn_analyse.pack(pady=5)
        btn_load_fasta = tk.Button(self.root,
                                   text="Load FASTA file",
                                   command=self.load_fasta_file)
        btn_load_fasta.pack(pady=5)
        label_result = tk.Label(self.root,
                                textvariable=self.result_text,
                                justify="left",
                                bg="white",
                                anchor="w",
                                width=80,
                                height=10,
                                relief="solid")
        label_result.pack(pady=10, padx=10)
        self.graph_frame = tk.Frame(self.root)
        self.graph_frame.pack(pady=10)

    def load_fasta_file(self):
        file_path = filedialog.askopenfilename(title="load FASTA file",
                                               filetypes=[("FASTA files",
                                                           "*.fasta *.fa"),
                                                          ("All files", "*.*")
                                                          ])
        if not file_path:
            return
        try:
            with open(file_path, 'r') as file:
                record = next(SeqIO.parse(file, "fasta"))
                self.dna_sourse = str(record.seq)
                self.entry_sequence.delete(0, tk.END)
                self.entry_sequence.insert(0, self.dna_sourse)
                self.analyze_sequence()
        except Exception as e:
            messagebox.showerror("Error", f"failed to upload the file: {e}")

    def analyze_sequence(self):
        self.dna_sourse = self.entry_sequence.get().upper()
        dna_analyzer = DNAAnalyzer(self.dna_sourse)
        if not dna_analyzer.validate_sequence():
            messagebox.showerror(
                "Error", "The sequence contains unacceptable characters")
            return
        counts = dna_analyzer.count_nucleotides()
        gc_content = dna_analyzer.gc_content()
        rna_sequence = dna_analyzer.transcribe()
        start_codons = dna_analyzer.find_codons(codons_type="start")
        stop_codons = dna_analyzer.find_codons(codons_type="stop")
        protein = dna_analyzer.translate_rna_to_protein(rna_sequence)
        self.result_text.set(
            f"Length of sequence: {len(self.dna_sourse)}\n"
            f"A: {counts['A']}, T: {counts['T']}, G: {counts['G']}, C: {counts['C']}\n"
            f"GC-content: {gc_content:.2f}%\n"
            f"Start-codons on positions: {start_codons}\n"
            f"Stop-codons on positions: {stop_codons}\n"
            f"RNA sequence: {rna_sequence}\n"
            f"Protein sequence: {protein}")
        self.plot_nucleotide_distribution(counts)

    def plot_nucleotide_distribution(self, counts):
        for widget in self.graph_frame.winfo_children():
            widget.destroy()
        labels = ['A', 'C', 'T', 'G']
        sizes = [counts['A'], counts['C'], counts['T'], counts['G']]
        color = ['#ff9999', '#66b3ff', '#99ff99', '#ffcc99']
        fig, ax = plt.subplots(figsize=(4, 4), dpi=100)
        ax.pie(sizes,
               labels=labels,
               autopct='%1.1f%%',
               startangle=90,
               colors=color)
        ax.axis('equal')
        canvas = FigureCanvasTkAgg(fig, master=self.graph_frame)
        canvas.draw()
        canvas.get_tk_widget().pack()


if __name__ == "__main__":
    root = tk.Tk()
    app = DNAApp(root)
    root.mainloop()
