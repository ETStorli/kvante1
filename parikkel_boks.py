from scipy.sparse import diags
import scipy.linalg as la
import numpy as np
import matplotlib.pyplot as plt
import scipy as sp
import scipy.constants as sc


l = 1e-9                      # Lengde på boksen
Dx = 1.0E-11                  # Steglengden
x = np.arange(0, l, Dx)       # x-kooridnatene
size = len(x)                 # Størrelsen på matrisa
V = np.zeros(size)            # Potensialet, som et array

# konstanter, definert ved hjelp av scipy - biblioteket
h = sc.h
hbar = sc.hbar
m = sc.m_e
c = sc.c

def analy(i):
    E = (i**2*sc.pi**2*hbar**2)/(2*m*l**2)
    return E/sc.eV


n = np.arange(1, (size+1), 1)
V = np.piecewise(x, [x < 0, x > l], [np.inf, np.inf, 0])
energi_funk = np.sqrt(2/l) - np.sin((n*np.pi*x)/l)
engeri_verdi = (n*np.pi*hbar)/(2*m*l)

