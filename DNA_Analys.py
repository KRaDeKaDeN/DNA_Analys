import tkinter as tk
from tkinter import filedialog, messagebox
from Bio import SeqIO
import matplotlib.pyplot as plt
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg


def analyse_dna(dna_sourse):
    dna_sourse = dna_sourse.upper()
    if not all(base in "ATGC" for base in dna_sourse):
        messagebox.showerror("Error, sequence contains invalid characters")
        return None
    count_a = dna_sourse.count("A")
    count_t = dna_sourse.count("T")
    count_g = dna_sourse.count("G")
    count_c = dna_sourse.count("C")

    gc_content = (count_g + count_c) / len(dna_sourse) * 100

    rna_sequence = dna_sourse.replace("T", "U")

    return {
        "length": len(dna_sourse),
        "count_a": count_a,
        "count_t": count_t,
        "count_g": count_g,
        "count_c": count_c,
        "GC_content": gc_content,
        "RNA_sequence": rna_sequence,
    }


def process_analysis(dna_sourse):
    results = analyse_dna(dna_sourse)
    if results is None:
        return
    result_text.set(
        f"The length of sequence is:{results['length']}\n"
        f"A= {results['count_a']}, T= {results['count_t']} G= {results['count_g']} C= {results['count_c']}\n"
        f"GC content= {results['GC_content']}\n"
        f"RNA sequence= {results['RNA_sequence']}"
    )

    plot_distr_result(results)


def load_Fasta_file():
    file_path = filedialog.askopenfilename(
        title="chuse FASTA file",
        filetypes=[("FASTA files*fasta *.fa"), ("All files", "*.*")],
    )
    if not file_path:
        return
    try:
        with open(file_path, "r") as file:
            record = next(SeqIO.parse(file, "fasta"))
            dna_sourse = str(record.seq)
            entry_sequence.delete(0, tk.END)
            entry_sequence.insert(0, dna_sourse)
            process_analysis(dna_sourse)
    except Exception as e:
        messagebox.showerror(f"Error failed to load file {e}")


def plot_distr_result(results):
    for widget in graph_frame.winfo_children():
        widget.destroy()
    labels = ["A", "T", "G", "C"]
    sizes = [
        results["count_a"],
        results["count_t"],
        results["count_g"],
        results["count_c"],
    ]
    colors = ["#ff9999", "#66b3ff", "#99ff99", "#ffcc99"]
    fig, ax = plt.subplots(figsize=(4, 4), dpi=100)
    ax.pie(sizes, labels=labels, autopct="%1.1f%%", startangle=90, colors=colors)
    ax.axis("equal")
    canvas = FigureCanvasTkAgg(fig, master=graph_frame)
    canvas.draw()
    canvas.get_tk_widget().pack()


def analyse_from_text():
    dna_sourse = entry_sequence.get().upper()
    process_analysis(dna_sourse)


window = tk.Tk()
window.title("DNA and RNA analys")
window.geometry("700x700")

label_instruction = tk.Label(window, text="Write DNA sequence or load FASTA file")
label_instruction.pack(pady=10)
entry_sequence = tk.Entry(window, width=80)
entry_sequence.pack(pady=10)

btn_analyse = tk.Button(window, text="Analyse", command=analyse_from_text)
btn_analyse.pack(pady=10)

btn_load_fasta = tk.Button(window, text="Load FASTA file")
btn_load_fasta.pack(pady=10)

result_text = tk.StringVar()
label_result = tk.Label(
    window,
    textvariable=result_text,
    justify="left",
    bg="white",
    anchor="w",
    width=60,
    height=10,
    relief="solid",
)
label_result.pack(pady=10, padx=10)

graph_frame = tk.Frame(window)
graph_frame.pack(pady=10)

window.mainloop()
