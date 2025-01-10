import tkinter as tk
from tkinter import messagebox

def analyse_dna():
    dna_sourse = entry_sequence.get().upper() 
    if not all(base in "ATGC" for base in dna_sourse):
        messagebox.showerror("Error, sequence contains invalid characters")
        return
    count_a = dna_sourse.count('A')
    count_t = dna_sourse.count('T')
    count_g = dna_sourse.count('G')
    count_c = dna_sourse.count('C')
    
    gc_content = ((count_g + count_c)) / len(dna_sourse) * 100
    
    rna_sequence = dna_sourse.replace('T', 'U')
    
    result_text.set(
            f"The length of sequence is: {len(dna_sourse)})\n"
            f"A: {count_a}, T: {count_t} G: {count_g} C: {count_c}\n"
            f"GC-content:{gc_content:.2f}%\n"
            f"RNA sequence: {rna_sequence}"
        )

window = tk.Tk()
window.title("DNA and RNA analys")
window.geometry("500x300")

label_instruction = tk.Label(window, text="Write DNA sequence:")
label_instruction.pack(pady=10)
entry_sequence = tk.Entry(window, width=50)
entry_sequence.pack(pady=10)

btn_analyse = tk.Button(window, text="Analyse", command=analyse_dna)
btn_analyse.pack(pady=10)

result_text = tk.StringVar()
label_result = tk.Label(window, textvariable=result_text, justify="left", bg="white", anchor='w', width=60, height=10, relief="solid")
label_result.pack(pady=10, padx=10)

window.mainloop()