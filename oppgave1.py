from scipy.sparse import diags
import scipy as sp
import scipy.constants as sc
import scipy.linalg as la
import numpy as np
import matplotlib.pyplot as plt

size = 3

l = 2*np.pi
x = np.linspace(0, l, size)
#h = sc.h
#hbar = sc.hbar
#m = sc.m_e
h = 1
hbar = 1
m = 1
Dx = l/len(x)
c = 1
V = np.identity(size)


H_dia = np.power(hbar, 2)/(m*np.power(Dx, 2)) + V[0][0]
H_subDia = -np.power(hbar, 2)/(2*m*np.power(Dx,2)) #−ℏ2/(2𝑚Δ𝑥2)

H = diags([H_subDia, H_dia, H_subDia], [-1, 0, 1], shape=(size, size)).toarray()
egen_verdi = la.eig(H)
print(egen_verdi[1])
