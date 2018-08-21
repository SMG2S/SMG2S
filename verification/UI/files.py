#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jul  9 16:22:11 2018

@author: qpetit
"""

import re
import sys
if sys.version_info[0] < 3:
    import tkFileDialog as filedialog
    import tkMessageBox as messagebox
else:
    from tkinter import messagebox
    from tkinter import filedialog
from decimal import Decimal
import numpy as np


from get_path import read_file_path
from graph import generate_graph

#Fonction qui permet la sauvegarde du graphique là ou l'utilisateur le veut
def save_file(b):
    try:
        fichier = filedialog.asksaveasfilename(filetypes=[("Portable Network Graphics", "*.png"),("Encapsulated Postcript", "*.eps"),("Portable Document Format", "*.pdf"),("Scalable Vector Graphics", "*.svg,*,svgz"),("Joint Photographic Experts Group", "*.jpeg,*.jpg"),("PGF code for LaTeX", "*.pgf"),("Raw RGBA bitmap", "*.raw,*.rgba"),("Postscript", "*.ps")])
    except FileNotFoundError:
        messagebox.showerror("error", "Problem with file path, the image has not been saved")
        return
    
    #Si l'utilisateur a sélecxtionné un chemin
    if fichier != "" :
        dimension = read_files()
        if dimension < 1 :
            return
        graphique = generate_graph(dimension, fichier, b)
        if graphique != True :
            messagebox.showerror("Error","Can not display the graphic")
    
#Fonction qui permet la lecture des fichiers .txt de données
def read_files():
    #On ouvre le fichier 1
    fichier1 = read_file_path(1)
    #On ouvre le fichier 2
    fichier2 = read_file_path(2)
    
    exp = r"[0-9]+(\.)?[0-9]*"
    # Booléen qui permet de stocker si on a déjà rencontré la ligne de la dimension.
    premiere_ligne = False
    dimension = 0
    if fichier1 == "" and fichier2 == "" :
        messagebox.showerror("Error","No files were opened")
        return -1
    else:
        if fichier1 != "":
            presence_fichier1 = True
            with open(fichier1, "r") as f:
                for line in f.readlines():
                    donnees = line.split()
                    
                    #On regarde si le premier caractère n'est pas un espace est un chiffre
                    if re.match(exp, donnees[0]) is not None:
                        if premiere_ligne == False:
                            premiere_ligne = True
                            dimension = 0
                            for x in line.split(" "):
                                if re.match(exp, x):
                                    dimension += 1
                            dimension = dimension - 1
                            #On créé les vecteurs qui stockera les informations a afficher
                            reel1 = []
                            imaginaire1 = []
                            
                        else:
                            reel1.append(Decimal(donnees[1]))
                            if dimension == 2:
                                imaginaire1.append(Decimal(donnees[2]))
                            else:
                                imaginaire1.append(0)
            f.close()    
            #On sauvegarde les données
            np.savetxt(r'data/matrix/r1.vec', reel1)
            np.savetxt(r'data/matrix/i1.vec', imaginaire1)
        else:      
            presence_fichier1 = False 
        
        
        premiere_ligne = False
        if fichier2 != "":
            presence_fichier2 = True
            with open(fichier2, "r") as f:
                for line in f.readlines():
                    donnees = line.split()
                    
                    #On regarde si le premier caractère n'est pas un espace est un chiffre
                    if re.match(exp, donnees[0]) is not None:
                        if premiere_ligne == False:
                            premiere_ligne = True
                            ancienne_dimension = dimension
                            dimension = 0
                            for x in line.split(" "):
                                if re.match(exp, x):
                                    dimension += 1
                            dimension = dimension - 1
                            #On créé les vecteurs qui stockera les informations a afficher
                            reel2 = []
                            imaginaire2 = []
                            
                        else:
                            reel2.append(Decimal(donnees[1]))
                            if dimension == 2:
                                imaginaire2.append(Decimal(donnees[2]))
                            else:
                                imaginaire2.append(0)
            f.close()
            #On sauvegarde les données
            np.savetxt(r'data/matrix/r2.vec', reel2)
            np.savetxt(r'data/matrix/i2.vec', imaginaire2)
        else:  
            presence_fichier2 = False
        
        #on va vérifier que les deux fichiers sont de la même taille
        if presence_fichier1 == presence_fichier2 == True :
            if ancienne_dimension != dimension :
                messagebox.showwarning("Error","The selected files are not the same size. Please check the files") 
                return -1
            else:
                return dimension
        else:
            return dimension


def open_file(number):
    
    #On charge les chemins déjà présent
    if number == 1 :
        fichier = read_file_path(1)
    else:
        fichier = read_file_path(2)
    #On enregistre ce que contient l'adresse du fichier précédente
    fichier_precedent = fichier
    
    #On demande à l'utilisateur de choisir son fichier texte
    fichier = filedialog.askopenfilename(filetypes=[('Text file','.txt')])
    
    #Si le fichier choisi est différent du fichier précédent alors le graphique n'est plus valable, on va donc pouvoir en refaire un et mettre à jour l'affichage.
    if fichier_precedent != fichier :
        file = open("data/files/file"+str(number)+".txt", "w")
        file.write(fichier)
    