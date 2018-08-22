#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jul  3 12:26:49 2018

@author: qpetit
"""


import sys
if sys.version_info[0] < 3:
    from Tkinter import * 
    import tkMessageBox as messagebox
else:
    from tkinter import *
    from tkinter import messagebox


def read_resolution():
    with open(r"sys/save.reso", "r") as f:
        for line in f.readlines():
            donnees = line.split()
            resolution = donnees[0]+" x "+donnees[1]
    f.close()
    return resolution

def save_setting_save(taille, fenetreParametre, choix_x, choix_y):
    fichier = open(r"sys/save.reso", "w")
    if taille.get() != "perso":
        fichier.write(taille.get())
    else:
        fichier.write(choix_x.get()+" "+choix_y.get())
    fichier.close()
    fenetreParametre.destroy()

def setting_save(fenetre):
    fenetreParametre=Toplevel(fenetre)
    fenetreParametre.title("Settings - SMG2S")
    taille = StringVar()
    var_abs = StringVar()
    var_ord = StringVar()
    titre_parametre = Label(fenetreParametre, text="selection of the saved chart resolution")
    choix_2 = Radiobutton(fenetreParametre, text="2 Mpx (1632 x 1224)", variable=taille, value="1632 1224")
    choix_10 = Radiobutton(fenetreParametre, text="10 Mpx (3648 x 2736)", variable=taille, value="3648 2736")
    choix_25 = Radiobutton(fenetreParametre, text="25 Mpx (5776 x 4336)", variable=taille, value="5776 4336")
    choix_48 = Radiobutton(fenetreParametre, text="48 Mpx (8000 x 6000)", variable=taille, value="8000 6000")
    choix_perso = Radiobutton(fenetreParametre, text="Custom : ", variable=taille, value="perso")
    choix_x = Label(fenetreParametre, text="W")
    choix_y = Label(fenetreParametre, text=" & H :")
    resolution = read_resolution()
    current = Label(fenetreParametre, text="current resolution : "+resolution)
    button_cancel = Button(fenetreParametre, text="Cancel", command=fenetreParametre.destroy)
    button_cancel.grid(row=7, column=1, columnspan = 2, sticky="w")
    button_save = Button(fenetreParametre, text="Save",command=lambda:save_setting_save(taille, fenetreParametre, var_abs, var_ord))
    button_save.grid(row=7, column=3, columnspan = 2, sticky="w")

    
    ligne_abs = Entry(fenetreParametre, textvariable=var_abs, width=10)
    
    ligne_ord = Entry(fenetreParametre, textvariable=var_ord, width=10)
    titre_parametre.grid(row=0, column=0, columnspan = 5, sticky="w")
    choix_2.grid(row=1, column=0, columnspan = 5, sticky="w")
    choix_10.grid(row=2, column=0, columnspan = 5, sticky="w")
    choix_25.grid(row=3, column=0, columnspan = 5, sticky="w")
    choix_48.grid(row=4, column=0, columnspan = 5, sticky="w")
    choix_perso.grid(row=5, column=0, columnspan = 1, sticky="w")
    choix_x.grid(row=5, column=1, columnspan = 1, sticky="w")
    ligne_abs.grid(row=5, column=2, columnspan = 1, sticky="w")
    choix_y.grid(row=5, column=3, columnspan = 1, sticky="w")
    ligne_ord.grid(row=5, column=4, columnspan = 1, sticky="w")
    current.grid(row=6, column=0, columnspan = 5, sticky="w")
    
    fenetreParametre.mainloop()
    
#------------------------------------------------------------------------------
if __name__ == '__main__':
    fenetre = Tk()
    setting_save(fenetre)
    fenetreParametre.mainloop()