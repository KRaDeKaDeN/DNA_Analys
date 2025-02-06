import tkinter as tk
from bio_analyzer.gui.main_window import BioApp


def start():
	root = tk.Tk()
	app = BioApp(root)  # noqa: F841
	root.mainloop()


if __name__ == "__main__":
	start()
