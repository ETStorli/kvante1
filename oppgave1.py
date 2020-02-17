from scipy.sparse import diags
import scipy.linalg as la
import numpy as np
import matplotlib.pyplot as plt
import scipy as sp
import scipy.constants as sc

size = 10000
l = 1e-9
dx = 1e-11
x = np.arange(0, l, dx)
number_of_funk = range(5)
h = sc.h
hbar = sc.hbar
m = sc.m_e
#h = 1
#hbar = 1
#m = 1
c = sc.c
#V = np.identity(size)
V = np.zeros(size)


H_dia = np.power(hbar, 2)/(m*np.power(dx, 2)) + V
H_subDia = [-np.power(hbar, 2)/(2*m*np.power(dx,2))]*(size-1)

#Bruker eigh istedet for eig. Vi kan bruke eigh siden H er symmetrisk, den er kjappere og sorterer egenverdiene
egenverdi, egenvektor = la.eigh_tridiagonal(H_dia, H_subDia)        #Gjør samme jobb som det som er kommentert ut ovenfor, gir samme verdier

egenvektor_abs_pow = np.power(np.abs(egenvektor), 2)
integral = np.sum(egenvektor_abs_pow)
normert = 1/la.norm(integral)
print(normert)

'''
plt.isinteractive()
plt.figure('Test',dpi=100, figsize=(13,5))
plt.subplot(121)
for i in number_of_funk:
    plt.plot(x, egenvektor[:, i])
plt.subplot(122)
for i in number_of_funk:
    plt.plot(x, [egenverdi[i]]*len(x))'''
#plt.show()
'''
for i in range(3, 6):
    plt.plot(x, egenvektor[:, i])
plt.show()

'''
print(egenverdi[0]/sc.eV)
print(egenverdi[1]/sc.eV)
print(egenverdi[2]/sc.eV)
print(egenverdi[3]/sc.eV)
