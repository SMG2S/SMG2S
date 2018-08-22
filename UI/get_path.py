#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jul  9 16:12:19 2018

@author: qpetit
"""

#Fonction qui retourne le chemin du fichier à ouvrir (Fichier1 et fichier2)
def read_file_path(number):
    path = "data/files/file"+str(number)+".txt"
    #On ouvre le fichier 1
    f = open(path, "r")
    fichier = f.read()
    f.close()

    if len(fichier) ==0:
        fichier = ""
    return fichier