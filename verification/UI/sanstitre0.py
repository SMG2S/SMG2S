#/usr/bin/python
# -*- coding: utf-8 -*-
 
import tkinter as tk
 
class TkinterApp:
 
    def __init__(self):
        # fenêtre principale de l'application :
        self.root = tk.Tk()
        self.root.title("main window")
        self.root.geometry('300x300+200+200')
 
        # bouton pour ajouter une nouvelle fenêtre
        self.bt = tk.Button(self.root, text="New window", command=self.newWindow)
        self.bt.grid(row=0, column=0)
 
    def newWindow(self):
        # créer une nouvelle fenêtre en cliquant sur le bouton
        self.newWindow = tk.Tk()
        self.newWindow.title("new window")
        self.newWindow.geometry('200x200+250+250')
        self.newWindow.mainloop()
 
# programme principal
if __name__ == '__main__':
 
     
    app = TkinterApp()
    app.root.mainloop()
    