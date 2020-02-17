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

from oppgave1 import numSchrodinger, plot_egenverdi, plot_egenvektor
w = 1
V_1 = 0.5*m*np.power(w, 2)*np.power(x, 2)
egenverdi, egenvektor = numSchrodinger(V_1)


plt.figure('Boks-potensial',dpi=100, figsize=(13,5))
plt.subplot(121)
plot_egenvektor(0, 3, egenvektor)

plt.subplot(122)
plot_egenverdi(4, egenverdi, True)

plt.show()

M = m*0.5


def Hermite(n,x):
    H = [1, 2*x]
    for i in range(2, n+1):
        H.append(2 * x * H[i - 1] - 2 * (i - 1) * H[i - 2])
    return H[n]


def analytisk_oscillator(n):
    psi = []
    L = 1 / 2
    X = np.linspace(-L * 1E-1, L * 1E-1, 101)
    for x in X:
        psi.append((1/(np.sqrt((2**n)*np.math.factorial(n))))*((M*w/(np.pi*sc.hbar))**(1/4)) *(np.exp(-(M*w*x**2)/(2*sc.hbar)))*(Hermite(n, np.sqrt(M*w/sc.hbar)*x)))
    return psi


def plot_egenvektor_liste(n):
    plt.figure()
    plt.title('Harmonisk oscillator, egenfunksjoner, analytisk')
    for i in range(n):
        plt.plot(x, analytisk_oscillator(i))
    plt.grid(True, linestyle='--')


def plot_egenverdi_analytisk(n):
    plt.title(r'Egenverdier $E_j$', fontsize=15)
    for i in range(n):
        plt.plot(x, hbar*w*(i+0.5)*len(x), label=r'$j = $%.i'%(i))
    plt.yticks((hbar*w*(0+0.5), hbar*w*(1+0.5), hbar*w*(2+0.5), hbar*w*(3+0.5)),
               ('%.2f eV'%(hbar*w*(0+0.5)), '%.2f eV'%(hbar*w*(1+0.5)),
                '%.2f eV'%(hbar*w*(2+0.5)), '%.2f eV'%(hbar*w*(3+0.5))))
    plt.xticks((), ())

plt.subplot(121)
plot_egenvektor_liste(4)

plt.subplot(122)
plot_egenverdi_analytisk(4)
plt.show()
