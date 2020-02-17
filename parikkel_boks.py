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
#h = sc.h
#hbar = sc.hbar
#m = sc.m_e
h = 1
hbar = 1
m = 1
Dx = l/len(x)
c = 1
V = np.identity(size)
