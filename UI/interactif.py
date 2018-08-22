#!/usr/bin/env python3
# -*- coding: utf-8 -*-
import sys
if sys.version_info[0] < 3:
    import Tkinter as Tk
else:
    import tkinter as Tk


import matplotlib
matplotlib.use("tkAgg")
from matplotlib.pyplot import scatter, legend, savefig, show, xlim, ylim

from numpy import loadtxt
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg, NavigationToolbar2TkAgg
# implement the default mpl key bindings
from matplotlib.backend_bases import key_press_handler

from matplotlib.figure import Figure


from get_path import read_file_path

#Fonction qui affiche le graphique demandé dans la mee fenêtre (dans la partie droite de la fenêtre)
def display_graph_ext(dimension, b=False):

    root = Tk.Tk()
    root.wm_title("SMG2S UI")
    
    
    graphique = True
    
    #On ouvre les fichiers contenant les chemins des deux fichiers selectionnes
    #On ouvre le fichier 1
    fichier1 = read_file_path(1)
    #On ouvre le fichier 2
    fichier2 = read_file_path(2)
    
    print(fichier1, " & ", fichier2)
    
    f = Figure(figsize=(8, 6), dpi=100)
    a = f.add_subplot(1, 1, 1)
    
    reel1= loadtxt(r'data/matrix/r1.vec')
    imaginaire1 = loadtxt(r'data/matrix/i1.vec')
    reel2 = loadtxt(r'data/matrix/r2.vec')
    imaginaire2 = loadtxt(r'data/matrix/i2.vec')
    
    
    if dimension == 2:
        if fichier1 != "":
            a.scatter(reel1, imaginaire1, c = 'black', marker = 'o', s = 200, label="Initial Eigenvalues")
            graphique = True
        if fichier2 != "":
            a.scatter(reel2, imaginaire2, c = 'red', marker = '+', s = 200, label="Computed Eigenvalues")
            graphique = True
        
    elif dimension == 1:
            
            if fichier1 != "":
                a.scatter(reel1, imaginaire1, c = 'black', marker = 'o', s = 200, label="Initial Eigenvalues")
                graphique = True
            if fichier2 != "":
                a.scatter(reel2, imaginaire2, c = 'red', marker = '+', s = 200, label="Computed Eigenvalues")
                graphique = True
    else: 
        Tk.messagebox.showerror("Error", "Impossible to generate the graph, the dimension is incorrect")
        graphique = False
    if graphique == True :
        if b == True:
            u = loadtxt(r'data/custom/axe.vec')
            a.set_xlim(u[0], u[1])
            a.set_ylim(u[2], u[3])
        a.legend(bbox_to_anchor=(0, 1.02, 1, .102), loc=3, ncol=2, borderaxespad=0)
        show()
    
    
    
    # a tk.DrawingArea
    canvas = FigureCanvasTkAgg(f, master=root)
    canvas.get_tk_widget().pack(side=Tk.TOP, fill=Tk.BOTH, expand=1)
    
    toolbar = NavigationToolbar2TkAgg(canvas, root)
    toolbar.update()
    canvas._tkcanvas.pack(side=Tk.TOP, fill=Tk.BOTH, expand=1)
    
    
    def on_key_event(event):
        print('you pressed %s' % event.key)
        key_press_handler(event, canvas, toolbar)
    
    canvas.mpl_connect('key_press_event', on_key_event)
    
    
    def _quit():
        root.quit()     # stops mainloop
        root.destroy()  # this is necessary on Windows to prevent
                        # Fatal Python Error: PyEval_RestoreThread: NULL tstate
    
    button = Tk.Button(master=root, text='Quit', command=_quit)
    button.pack(side=Tk.BOTTOM)
    
    Tk.mainloop()
# If you put root.destroy() here, it will cause an error if
# the window is closed with the window manager.
    
    

#Fonction qui affiche le graphique demandé dans la mee fenêtre (dans la partie droite de la fenêtre)
def display_graph_int(dimension, root, canvas, b=False):
    
    for w in root.winfo_children():
        w.destroy()
    graphique = True
    
    #On ouvre les fichiers contenant les chemins des deux fichiers selectionnes
    #On ouvre le fichier 1
    fichier1 = read_file_path(1)
    #On ouvre le fichier 2
    fichier2 = read_file_path(2)
    
    print(fichier1, " & ", fichier2)
    
    f = Figure(figsize=(8, 6), dpi=100)
    a = f.add_subplot(1, 1, 1)
    
    reel1= loadtxt(r'data/matrix/r1.vec')
    imaginaire1 = loadtxt(r'data/matrix/i1.vec')
    reel2 = loadtxt(r'data/matrix/r2.vec')
    imaginaire2 = loadtxt(r'data/matrix/i2.vec')
    
    
    if dimension == 2:
        if fichier1 != "":
            a.scatter(reel1, imaginaire1, c = 'black', marker = 'o', s = 200, label="Initial Eigenvalues")
            graphique = True
        if fichier2 != "":
            a.scatter(reel2, imaginaire2, c = 'red', marker = '+', s = 200, label="Computed Eigenvalues")
            graphique = True
        
    elif dimension == 1:
            
            if fichier1 != "":
                a.scatter(reel1, imaginaire1, c = 'black', marker = 'o', s = 200, label="Initial Eigenvalues")
                graphique = True
            if fichier2 != "":
                a.scatter(reel2, imaginaire2, c = 'red', marker = '+', s = 200, label="Computed Eigenvalues")
                graphique = True
    else: 
        Tk.messagebox.showerror("Error", "Impossible to generate the graph, the dimension is incorrect")
        graphique = False
    if graphique == True :
        if b == True:
            u = loadtxt(r'data/custom/axe.vec')
            a.set_xlim(u[0], u[1])
            a.set_ylim(u[2], u[3])
        a.legend(bbox_to_anchor=(0, 1.02, 1, .102), loc=3, ncol=2, borderaxespad=0)
        show()
    
    
    
    # a tk.DrawingArea
    canvas = FigureCanvasTkAgg(f, master=root)
    canvas.get_tk_widget().pack(side=Tk.TOP, fill=Tk.BOTH, expand=1)
    
    toolbar = NavigationToolbar2TkAgg(canvas, root)
    toolbar.update()
    canvas._tkcanvas.pack(side=Tk.TOP, fill=Tk.BOTH, expand=1)
    
    
    def on_key_event(event):
        print('you pressed %s' % event.key)
        key_press_handler(event, canvas, toolbar)
    
    canvas.mpl_connect('key_press_event', on_key_event)
    
    
    
    def _quit():
#        canvas.quit()     # stops mainloop
#        canvas.destroy()  # this is necessary on Windows to prevent
                        # Fatal Python Error: PyEval_RestoreThread: NULL tstate
        for w in root.winfo_children():
            w.destroy()
    
    button.pack(side=Tk.BOTTOM)
    
    Tk.mainloop()