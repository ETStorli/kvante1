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


def numSchrodinger(V, size = size):         #V må være en array
    H_dia = np.power(hbar, 2) / (m * np.power(Dx, 2)) + V
    H_subDia = np.array([-np.power(hbar, 2) / (2 * m * np.power(Dx, 2))] * (size - 1))
#    H = diags([H_subDia, H_dia, H_subDia], [-1, 0, 1], shape=(size, size)).toarray() #Trengs kanskje ikke
    egenverdi, egenvektor = la.eigh_tridiagonal(H_dia, H_subDia)    #Gir oss egenverdiene og egenvektorene til matrisen H
    egenverdi *= 1 / sc.eV           #Gjør om til elektronvolt

    for i in range(len(egenvektor)):
        egenvektor_abs_pow = np.power(np.abs(egenvektor[i]), 2)
        integral = np.sum(egenvektor_abs_pow)             # "Integralet" av absoluttverdien til energiegenfunksjonene i andre
        normert_konst = np.sqrt(1 / la.norm(integral))    # Normeringskonstanten
        egenvektor[i] *= normert_konst

    return egenverdi, egenvektor


def plot_egenverdi(n, egenverdi, plotPotensial):
    plt.title(r'Egenverdier $E_j$', fontsize=15)
    if plotPotensial:           #Plotter potensialet
        plt.plot(x, V/1.60E-19, label=r'$V(x)$')
    for i in range(n):
        plt.plot(x, [egenverdi[i]]*len(x), label=r'$j = $%.i'%(i))
    plt.yticks((egenverdi[0], egenverdi[1], egenverdi[2], egenverdi[3]),
               ('%.2f eV'%(egenverdi[0]), '%.2f eV'%(egenverdi[1]),
                '%.2f eV'%(egenverdi[2]), '%.2f eV'%(egenverdi[3])), fontsize=10)
    plt.xticks((), ())


def plot_egenvektor(n1, n2, egenvektor):
    plt.title(r'Energiegenfunksjonene $\vec{\psi}_{(j)}$', fontsize=15)
    for i in range(n1, n2 + 1):
        plt.plot(x, egenvektor[:, i], label=r'$j = $%.i' % (i + 1))
    plt.xlabel('$x$ (nm)')
    plt.ylabel('$\psi(x)$')
    plt.grid(True)
    plt.legend(loc='lower left')


egenverdi, egenvektor = numSchrodinger(V)
plt.figure('Boks-potensial',dpi=100, figsize=(13,5))
plt.subplot(121)
plot_egenvektor(0, 3, egenvektor)

plt.subplot(122)
plot_egenverdi(4, egenverdi, False)        #Med argument_2 = False/True --> plotter ikke eller plotter potensialet

plt.show()

print(egenverdi[0], egenverdi[1], egenverdi[2])
