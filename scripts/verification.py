'''
MIT License
Copyright (c) 2019 Xinzhe WU @ Maison de la Simulation, France
Copyright (c) 2019-2022, Xinzhe Wu @ Simulation and Data Laboratory Quantum 
                                     Materials,  Forschungszentrum Juelich GmbH.
                                     
Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:
The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.
THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.
'''

import scipy.io
from tabulate import tabulate
import matplotlib.pyplot as plt
import numpy as np
from scipy.sparse.linalg import eigs, spsolve
import scipy.sparse as sp
import pandas as pd
import argparse
from matplotlib.patches import Rectangle
import matplotlib.colorbar as cbar
import matplotlib.cm as cm
import matplotlib.collections as collections

#load sparse matrix from MatrixMarket format
#output is in COO format
def loadMatrix(filename):
    M = scipy.io.mmread(filename)
    return M 

#load vector from files following the format of SMG2S
def loadVector(filename):
    df=pd.read_csv(filename, comment="%", delim_whitespace=True)
    n = df.shape[1]

    if n == 3:
        return df[df.columns[1]].to_numpy() + 1j * df[df.columns[2]].to_numpy()
    elif n == 2:
        return df[df.columns[1]].to_numpy()
    else:
        raise ValueError('Oops! The given vector file is not in good format')

#shifted inverse iteration method to approach eigenvalue closest to mu
def shiftInverse(A, mu, tol=0.001, max_iter = 100):
    m, n = A.shape

    x = np.random.rand(m)

    I = sp.identity(m)

    idx = 0

    eval = 0
    eval_old = 0

    while idx < max_iter:
        Ashift = A - (mu*I)
        y = spsolve(Ashift, x) 
        x = y / np.linalg.norm(y)
        c = y.dot(x) / x.dot(x)
        eval = 1 / c + mu
        if np.abs(eval - eval_old) < tol * np.abs(eval_old):
            break
        eval_old = eval
        idx = idx + 1

    return eval, idx

def spy_coo(M, ax, type="pattern"):
    if not isinstance(M, sp.coo_matrix):
        M = sp.coo_matrix(M)

    verts = [((x-0.5, y-0.5), (x-0.5, y+0.5), (x+0.5, y+0.5), (x+0.5,
y-0.5)) for (x,y) in zip(M.col, M.row)]

    c = collections.PolyCollection(verts)
    if type == "heatmap":
        c.set_array(np.absolute(M.data) )
    c.set_cmap(cm.Wistia)

    ax.add_collection(c)
    ax.set_xlim(-1,8)
    ax.set_ylim(-1,8)

    ax.set_xlim(-0.5, M.shape[1]-0.5)
    ax.set_ylim(-0.5, M.shape[0]-0.5)
    ax.invert_yaxis()
    ax.set_aspect(float(M.shape[0])/float(M.shape[1]))
    ax.set_xticks([])
    ax.set_yticks([])
    ax.set_title("sparsity pattern")


    return ax

def plot_spectrum(input, estimate, ax):

    x_in = []
    y_in = []
    x_e = []
    y_e = []

    for v in input:
        x_in.append(v.real)
        if np.iscomplex(v):
            y_in.append(v.imag)
        else:
            y_in.append(0.0)

    for v in estimate:
        x_e.append(v.real)
        if np.iscomplex(v):
            y_e.append(v.imag)
        else:
            y_e.append(0.0)

    ax.scatter(x_in, y_in, c='black')
    ax.scatter(x_e, y_e, marker="+",c='r')
    ax.set_ylabel('Imaginary')
    ax.set_xlabel('Real')
    asp = np.diff(ax2.get_xlim())[0] / np.diff(ax2.get_ylim())[0]
    ax.set_aspect(asp)
    ax.legend(['Given spectrum', 'Computed eigenvalues'])
    ax.set_title("spectrum")

    return ax

def approx_spectrum(M, spec, offset):
    M_csr = M.tocsr()
    eigenvalues=[]
    for i in range(M.shape[0]):
        eval, _=shiftInverse(M_csr, (1+offset) * spec[i])
        eigenvalues.append(eval)

    return eigenvalues    

if __name__ == '__main__':

    parser = argparse.ArgumentParser(description='verification of matrices generated matrices to keep given spectra')

    parser.add_argument("--offset", default=1e-10)
    parser.add_argument("--matpath", default="data/testmatrix_cmplx.mtx")
    parser.add_argument("--specpath", default="data/given_spectrum_cmplx.txt")
    

    value = parser.parse_args()

    M=loadMatrix(value.matpath)
    spec = loadVector(value.specpath)
    eigenvalues = approx_spectrum(M, spec, value.offset)

    fig = plt.figure()
    ax = fig.add_subplot(121)
    ax = spy_coo(M, ax, type="heatmap")
    ax2 = fig.add_subplot(122)
    ax2 = plot_spectrum(spec, eigenvalues, ax2)

    plt.tight_layout()
    
    plt.show()


