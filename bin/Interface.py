from tkinter import *
from tkinter import filedialog
import os

import Prepare
import BlastTools
import Verify
import Iteration


window = Tk()
window.title("Recombination Analysis")
window.geometry('1050x300')


Label(window, text='Prepare', font=('Arial Bold', 30)).grid(columnspan=9, row=0, sticky=W)

# build blast database
Label(window, text='Build database').grid(column=0, row=1, sticky=W)
Label(window, text='in:').grid(column=1, row=1)
db_in_path = StringVar()
def select_db():
    filepath = filedialog.askopenfilename()
    db_in_path.set(filepath)
Entry(window, textvariable=db_in_path).grid(column=2, row=1)
Button(window, text="...", command=select_db).grid(column=3, row=1)

Label(window, text='out: ').grid(column=4, row=1)
db_out_path = StringVar()
def select_db2():
    filepath = filedialog.askdirectory()
    db_out_path.set(filepath)
Entry(window, textvariable=db_out_path).grid(column=5, row=1)
Button(window, text="...", command=select_db2).grid(column=6, row=1)

Label(window, text='database name: ').grid(column=7, row=1)
db_name = Entry(window, width=10)
db_name.grid(column=8, row=1)
def build_db():
    Prepare.build_db(db_in_path.get(), os.path.join(db_out_path.get(), db_name.get()), input_type='fasta')
Button(window, text="Build database", command=build_db).grid(column=9, row=1)

# convert_to_binary
Label(window, text='Convert binary file').grid(column=0, row=2, sticky=W)
Label(window, text='in:').grid(column=1, row=2)
idfile_in_path = StringVar()
def select_idfile():
    filepath = filedialog.askopenfilename()
    idfile_in_path.set(filepath)
Entry(window, textvariable=idfile_in_path).grid(column=2, row=2)
Button(window, text="...", command=select_idfile).grid(column=3, row=2)

Label(window, text='out: ').grid(column=4, row=2)
idfile_out_path = StringVar()
def select_idfile2():
    filepath = filedialog.asksaveasfilename(defaultextension='.acc')
    idfile_out_path.set(filepath)
Entry(window, textvariable=idfile_out_path).grid(column=5, row=2)
Button(window, text="...", command=select_idfile2).grid(column=6, row=2)

def convert_to_binary():
    Prepare.convert_to_binary(idfile_in_path.get(), idfile_out_path.get())
Button(window, text="Convert file", command=convert_to_binary).grid(column=9, row=2)


Label(window, text='Main', font=('Arial Bold', 30)).grid(columnspan=9, row=4, sticky=W)

# query
Label(window, text='Query: ').grid(column=0, row=5)
query_path = StringVar()
def select_query():
    filepath = filedialog.askopenfilename()
    query_path.set(filepath)
Entry(window, textvariable=query_path).grid(column=1, row=5)
Button(window, text="...", command=select_query).grid(column=2, row=5)

# backbone
Label(window, text='Backbone: ').grid(column=3, row=5)
backbone_path = StringVar()
def select_backbone():
    filepath = filedialog.askopenfilename()
    backbone_path.set(filepath)
Entry(window, textvariable=backbone_path).grid(column=4, row=5)
Button(window, text="...", command=select_backbone).grid(column=5, row=5)

# negative_seqidlist
Label(window, text='negative_seqidlist: ').grid(column=6, row=5)
negative_path = StringVar()
def select_negative():
    filepath = filedialog.askopenfilename()
    negative_path.set(filepath)
Entry(window, textvariable=negative_path).grid(column=7, row=5)
Button(window, text="...", command=select_negative).grid(column=8, row=5)

# database
Label(window, text='Database: ').grid(column=0, row=6)
text_db = Entry(window, width=20)
text_db.grid(column=1, row=6)
Label(window, text='Window size: ').grid(column=3, row=6)
text_window_size = Entry(window, width=20)
text_window_size.grid(column=4, row=6)
Label(window, text='Num cpus: ').grid(column=6, row=6)
text_num_cpus = Entry(window, width=20)
text_num_cpus.grid(column=7, row=6)

# output
Label(window, text='Out Log: ').grid(column=0, row=7)
log_out_path = StringVar()
def select_log():
    filepath = filedialog.asksaveasfilename(defaultextension='.acc')
    log_out_path.set(filepath)
Entry(window, textvariable=log_out_path).grid(column=1, row=7)
Button(window, text="...", command=select_log).grid(column=2, row=7)

def run_main():
    query = query_path.get()
    backbone = backbone_path.get()
    negative_seqidlist = negative_path.get()
    database = os.path.join(db_out_path.get(), text_db.get())
    window_size = int(text_window_size.get())
    num_cpus = int(text_num_cpus.get())
    log_path = log_out_path.get()
    Iteration.main(log_path, query, negative_seqidlist, database, num_cpus, window_size, _backbone=backbone)
    
Button(window, text="RUN", command=run_main).grid(column=9, row=9)

window.mainloop()