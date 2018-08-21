#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jul 10 09:49:42 2018

@author: qpetit
"""

import matplotlib.pyplot as plt
import numpy as np

import sys
if sys.version_info[0] < 3:
    import tkMessageBox as messagebox
else:
    from tkinter import messagebox

from get_path import read_file_path

def read_resolution():
    with open(r"sys/save.reso", "r") as f:
        for line in f.readlines():
            donnees = line.split()
            w = int(donnees[0])/300
            h = int(donnees[1])/300
    f.close()
    return w, h

def generate_graph(dimension, fichier, b=False):
    
    graphique = True
    
    #On ouvre les fichiers contenant les chemins des deux fichiers selectionné
    #On ouvre le fichier 1
    fichier1 = read_file_path(1)
    #On ouvre le fichier 2
    fichier2 = read_file_path(2)
    
    #On importe les vecteurs de valeurs
    reel1= np.loadtxt(r'data/matrix/r1.vec')
    imaginaire1 = np.loadtxt(r'data/matrix/i1.vec')
    reel2 = np.loadtxt(r'data/matrix/r2.vec')
    imaginaire2 = np.loadtxt(r'data/matrix/i2.vec')
    w, h =read_resolution()
    axe = plt.figure(1, figsize=(w, h))
    a = plt.subplot(1, 1, 1)
    
    
    if dimension == 2:
        if fichier1 != "":
            plt.scatter(reel1, imaginaire1, c = 'black', marker = 'o', s = 200, label="Original solution")
            graphique = True
        if fichier2 != "":
            plt.scatter(reel2, imaginaire2, c = 'red', marker = '+', s = 200, label="Final solution")
            graphique = True
        
    elif dimension == 1:
            
            if fichier1 != "":
                plt.scatter(reel1, imaginaire1, c = 'black', marker = 'o', s = 200, label="Original solution")
                graphique = True
            if fichier2 != "":
                plt.scatter(reel2, imaginaire2, c = 'red', marker = '+', s = 200, label="Final solution")
                graphique = True
    else: 
        messagebox.showerror("Error", "Impossible to generate the graph, the dimension is incorrect")
        graphique = False
        
    #Si tout est ok pour l'enregistrement, on savegarde le graphique
    if graphique == True :
        if b == True:
            u = np.loadtxt(r'data/custom/axe.vec')
            a.set_xlim(u[0], u[1])
            a.set_ylim(u[2], u[3])
        plt.legend(loc='best')
        plt.savefig(fichier, dpi=300)
        plt.show()
    return graphique