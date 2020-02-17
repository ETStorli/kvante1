from scipy.sparse import diags
import scipy.linalg as la
import numpy as np
import matplotlib.pyplot as plt
import scipy as sp
import scipy.constants as sc

size = 100
l = 2*np.pi
x = np.linspace(0, l, size)
number_of_funk = range(3)
h = sc.h
hbar = sc.hbar
m = sc.m_e
#h = 1
#hbar = 1
#m = 1
Dx = l/len(x)
c = sc.c


n = np.arange(1, (size+1), 1)
V = np.piecewise(x, [x < 0, x > l], [np.inf, np.inf, 0])
energi_funk = np.sqrt(2/l) - np.sin((n*np.pi*x)/l)
engeri_verdi = (n*np.pi*hbar)/(2*m*l)

