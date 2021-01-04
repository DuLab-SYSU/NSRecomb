from tkinter import *
from tkinter import filedialog
import os

window = Tk()
window.title("Recombination Analysis")
window.geometry('350x200')

lbl = Label(window, text='Query: ')
lbl.grid(column=0, row=0)

def clicked():
    file = filedialog.askopenfilenames(initialdir=os.path.dirname(__file__))
    lbl.configure(text=file.name)
btn = Button(window, text="...", command=clicked)
btn.grid(column=1, row=0)




window.mainloop()