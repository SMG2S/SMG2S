#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jul  2 15:42:00 2018

@author: qpetit
"""

import PIL.Image
import PIL.ImageTk
import numpy as np
from tkinter import *
from tkinter import messagebox
from matplotlib.figure import Figure
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg, NavigationToolbar2TkAgg


#Importation des fichiers .py présent dans le dossier 
from get_path import read_file_path
from files import read_files, open_file, save_file
from reset import reset_window
from settings import setting_save
from interactif import display_graph_ext, display_graph_int


def alert():
    messagebox.showinfo("Say Hello", "Hello World")
    
def save_file_launcher():
    global var_custom
    print(var_custom)
    if var_custom.get() == 0:
        save_file(False)
    else:
        b = stockage_custom_axe()
        save_file(b)

def mise_a_jour_interface():
    fichier1 = read_file_path(1)
    fichier2 = read_file_path(2)
   
    
    #Nettoyage du canves affichant les adresses
#    can.delete(ALL)
    
    #Ici, on regarde si le chemin du fichier est vide, si c'est le cas, on affiche un message d'attention, sinon on affiche le chemin du fichier.
    if fichier1 != "":
        button_file1.configure(fg = "green")
        button_display.configure(state="normal")
        button_save.configure(state="normal")
        button_display_windows.configure(state="normal")
        
        menubar.entryconfigure(2, state=NORMAL)
        menubar.entryconfigure(3, state=NORMAL)
    else:
        button_file1.configure(fg = "#F3421C")
    if fichier2 != "":
        
        button_file2.configure(fg = "green")
        button_display.configure(state="normal")
        button_save.configure(state="normal")
        button_display_windows.configure(state="normal")
        
        menubar.entryconfigure(2, state=NORMAL)
        menubar.entryconfigure(3, state=NORMAL)
    else: 
        button_file2.configure(fg = "#F3421C")
    
    if fichier1 == "" and fichier2 == "":
        button_display.configure(state="disabled")
        button_save.configure(state="disabled")
        button_display_windows.configure(state="disabled")
        menubar.entryconfigure(2, state=DISABLED)
        menubar.entryconfigure(3, state=DISABLED)
    fenetre.update_idletasks()


    
    
def affichage_fichier():
    
    
    dimension = read_files()
    if dimension < 1 :
        return
    else:
        global var_custom
        print(var_custom)
        if var_custom.get() == 0:
            display_graph_int(dimension, dynamic_screen, canvas, False)
        else:
            b = stockage_custom_axe()
            display_graph_int(dimension, dynamic_screen, canvas, b)
    
#        graph_origin = PIL.Image.open('data/graphic_brut.png')
#        graph_reduit = graph_origin.resize((h,l))        
#        img = PIL.ImageTk.PhotoImage(graph_reduit)
#        dic['image']= img
#        item = can_plot.create_image(h/2, l/2, image=img)
#        mon_image=PhotoImage(file=r"smiley.gif")
#        dic['image']= mon_image
#        img=can_plot.create_image(250,250,image=mon_image)
    
    fenetre.update_idletasks()
    fenetre.update()
        
        
    
def ouvrir_fichier1():
    open_file(1)
    mise_a_jour_interface()
    
    
    
    
def ouvrir_fichier2():
    open_file(2)
    mise_a_jour_interface()





def reset_window_launcher():
    global dic
    reset_window(dynamic_screen, dic, input_xmin, input_xmax, input_ymin, input_ymax, check_custom, button_display, button_save, button_file1, button_file2, menubar, button_display_windows)
    
def about():
    messagebox.showinfo("About...", "Implementation of a user interface for the SMG2S project: https://github.com/brunowu/SMG2S")
    return

def set_canvas_launcher():
    setting_save(fenetre)

def display_graph_launcher():
    dimension = read_files()
    if dimension < 1 :
        return
    else:
        global var_custom
        print(var_custom)
        if var_custom.get() == 0:
            display_graph_ext(dimension, False)
        else:
            b = stockage_custom_axe()
            if b == True:
                display_graph_ext(dimension, True)
            else: 
                messagebox.showerror("Error", "The coordinates entered are not valid")
                display_graph_ext(dimension, False)
   
def stockage_custom_axe():
    global input_xmin
    global input_xmax
    global input_ymin
    global input_ymax
    print(input_xmin.get())
    print(input_xmax.get())
    print(input_ymin.get())
    print(input_ymax.get())
    ####################ERREUR HERE#######################
    
    try:
        if float(input_xmin.get()) < float(input_xmax.get()) and float(input_ymin.get()) < float(input_ymax.get()):
            u = [float(input_xmin.get()), float(input_xmax.get()), float(input_ymin.get()), float(input_ymax.get())]
            np.savetxt(r'data/custom/axe.vec', u)
            print("on est là")
            return True
        else:
            messagebox.showerror("Error", "Interval problem. Please check the values ​​entered in the entries")
            return False
    except ValueError:
        messagebox.showerror("Error", "Please check the values ​​entered in the entries, there is a value error")
        
        
    
        
        
def custom():
    global var_custom
    print(var_custom)
    if var_custom.get() == 0:
        input_xmin.configure(state="disabled")
        input_xmax.configure(state="disabled")
        input_ymin.configure(state="disabled")
        input_ymax.configure(state="disabled")
    else:
        input_xmin.configure(state="normal")
        input_xmax.configure(state="normal")
        input_ymin.configure(state="normal")
        input_ymax.configure(state="normal")

#On ouvre la fenetre graphique de TKinter
fenetre = Tk()
fenetre.title("SMG2S - Home Page")
fenetre.configure(bg = "white")

#L = fenetre.winfo_screenwidth()
#H = fenetre.winfo_screenheight()
#
#
#print("largeur = ", L, "Hauteur = ", H)
#fenetre.geometry("%dx%d%+d%+d" % (L,H,0,0))

dimension = 0
dic={}





menubar = Menu(fenetre)
#On créé le menu déroulant 1 : Fichier
menu1 = Menu(menubar, tearoff=0)
menu1.add_command(label="Open original file", command=ouvrir_fichier1)
menu1.add_command(label="Open final file", command=ouvrir_fichier2)
menu1.add_separator()
menu1.add_command(label="Run", command=affichage_fichier)
menu1.add_separator()
menu1.add_command(label="Reset", command=reset_window_launcher)
menu1.add_command(label="Exit", command=fenetre.quit)
menubar.add_cascade(label="File", menu=menu1)

menu2 = Menu(menubar, tearoff=0)
menu2.add_command(label="Open in a new window", command=display_graph_launcher)
menubar.add_cascade(label="Display", menu=menu2)

menu3 = Menu(menubar, tearoff=0)
menu3.add_command(label="Save", command=save_file_launcher)
menubar.add_cascade(label="Save", menu=menu3)

#On créé le menu déroulant 3 : paramètres
menu4 = Menu(menubar, tearoff=0)
menu4.add_command(label="Canvas", command=set_canvas_launcher)
menubar.add_cascade(label="Settings", menu=menu4)

#On créé le menu déroulant 4 : Aide
menu5 = Menu(menubar, tearoff=0)
menu5.add_command(label="About...", command=about)
menubar.add_cascade(label="Help", menu=menu5)

fenetre.config(menu=menubar)




left_tab = Frame(fenetre, bg="white")
left_tab.grid(row=0, column=0, sticky="N")

title_files = Label(left_tab, text="Files", bg="white", font="bold 18")
title_files.grid(row=0, column=0, columnspan = 4, sticky="w")

title_select_files = Label(left_tab, text="Select files for display your graph", bg="white", bd = "5", font="bold")
title_select_files.grid(row=1, column=0, columnspan = 4, sticky="w")

text_file1 = Label(left_tab, text="Original file", bg="white")
text_file1.grid(row=2, column=0, columnspan = 2, sticky="")

button_file1 = Button(left_tab, text="Select",command=ouvrir_fichier1, fg="#F3421C")
button_file1.grid(row=3, column=0, columnspan = 2, sticky="")

text_file2 = Label(left_tab, text="Final file", bg="white")
text_file2.grid(row=2, column=2, columnspan = 2, sticky="")

button_file2 = Button(left_tab, text="Select",command=ouvrir_fichier2, fg="#F3421C")
button_file2.grid(row=3, column=2, columnspan = 2, sticky="")


title_display = Label(left_tab, text="Display & Save", bg="white", font="bold 18")
title_display.grid(row=4, column=0, columnspan = 4, sticky="w", pady=(15,0))

title_display_auto = Label(left_tab, text="Scales are automatically managed.", bg="white", bd = "5")
title_display_auto.grid(row=5, column=0, columnspan = 4, sticky="w")

var_custom = IntVar()
check_custom = Checkbutton(left_tab, text="Make with custom scale", variable=var_custom, bg = "white", command=custom)
check_custom.grid(row=6, column=0, columnspan = 4, sticky="w")


text_xmin = Label(left_tab, text="xmin : ", bg="white")
text_xmin.grid(row = 7, column = 0, sticky="")

input_xmin = Entry(left_tab, bg = "white", width="5")
input_xmin.grid(row = 7, column = 1, sticky = "")

text_xmax = Label(left_tab, text="xmax : ", bg="white")
text_xmax.grid(row = 7, column = 2, sticky="")

input_xmax = Entry(left_tab, bg = "white", width="5")
input_xmax.grid(row = 7, column = 3, sticky = "")


text_ymin = Label(left_tab, text="ymin : ", bg="white")
text_ymin.grid(row = 8, column = 0, sticky="")

input_ymin = Entry(left_tab, bg = "white", width="5")
input_ymin.grid(row = 8, column = 1, sticky = "")

text_ymax = Label(left_tab, text="ymax : ", bg="white")
text_ymax.grid(row = 8, column = 2, sticky="")

input_ymax = Entry(left_tab, bg = "white", width="5")
input_ymax.grid(row = 8, column = 3, sticky = "")

text_ymax = Label(left_tab, text=" ", bg="white", font = "1")
text_ymax.grid(row = 9, column = 0, sticky="")

button_display = Button(left_tab, text="Display",command=affichage_fichier)
button_display.grid(row=10, column=0, columnspan = 2, sticky="")



button_display_windows = Button(left_tab, text="New window",command=display_graph_launcher)
button_display_windows.grid(row=10, column=2, columnspan = 2, sticky="")

text_ymax = Label(left_tab, text=" ", bg="white", font = "1")
text_ymax.grid(row = 11, column = 0, sticky="")

title_display_auto = Label(left_tab, text="To save the chart with the quality selected", bg="white", bd = "5")
title_display_auto.grid(row=12, column=0, columnspan = 4, sticky="w")

button_save = Button(left_tab, text="Save",command=save_file_launcher)
button_save.grid(row=13, column=0, columnspan = 4, sticky="")



#Affichage de la partie droite de la fenêtre: l'affichage graphique dynamique 
dynamic_screen = Frame(fenetre, bg="white")
dynamic_screen.grid(row=0, column=1, sticky="NSEW")


fenetre.grid_columnconfigure(1, weight=1)
fenetre.grid_rowconfigure(0, weight=1)

canvas = Canvas(dynamic_screen,width=800, height=600, bg='white')
canvas.grid(sticky="NSEW")

reset_window_launcher()
fenetre.update_idletasks()
fenetre.update()

#ouverture_menu(fenetre)

# On démarre la boucle Tkinter qui s'interompt quand on ferme la fenêtre
fenetre.mainloop()

