from scipy.sparse import diags
import scipy.linalg as la
import numpy as np
import matplotlib.pyplot as plt
import scipy as sp
import scipy.constants as sc

size = 4
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


H_dia = [np.power(hbar, 2)/(m*np.power(Dx, 2)) + V[0][0]]*size
H_subDia = [-np.power(hbar, 2)/(2*m*np.power(Dx,2))]*(size-1) #−ℏ2/(2𝑚Δ𝑥2)
'''
H = diags([H_subDia, H_dia, H_subDia], [-1, 0, 1], shape=(size, size)).toarray()
#Bruker eigh istedet for eig. Vi kan bruke eigh siden H er symmetrisk, den er kjappere og sorterer egenverdiene
egen_verdi = la.eigh(H)
'''
egenverdi, egenvektor = la.eigh_tridiagonal(H_dia, H_subDia)        #Gjør samme jobb som det som er kommentert ut ovenfor, gir samme verdier

print(egen_verdi)
print(egenverdi)
print(egenvektor)

