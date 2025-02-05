import tkinter as tk
from .gui.main_window import BioApp

if __name__ == "__main__":
    root = tk.Tk()
    app = BioApp(root)
    root.mainloop()