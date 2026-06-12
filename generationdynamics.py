# -*- coding: utf-8 -*-
"""
Created on Wed Apr  1 21:51:02 2026

@author: Bente
"""


import numpy as np
import matplotlib.pyplot as plt

h_S = 0.5
c_sup = 0.15
h_D = 0.4

k = 0.6
max_generations = 5000

def next_gen(x00, x01, x10, x11, h_D, h_S, t, c_trait, c_sup):
    f_x00 = (
        x00**2
        + x00*x01
        + (1-h_D*t)*(1-h_D*c_trait)*x00*x10
        + 0.5*(1-h_D*(1-h_S)*t)*(1-h_D*(1-h_S)*c_trait - h_S*c_sup)*x00*x11
        + 0.5*(1-h_D*(1-h_S)*t)*(1-h_D*(1-h_S)*c_trait - h_S*c_sup)*x01*x10
    )

    f_x01 = (
        x00*x01
        + 0.5*(1-h_D*(1-h_S)*t)*(1-h_D*(1-h_S)*c_trait - h_S*c_sup)*x00*x11
        + x01**2
        + 0.5*(1-h_D*(1-h_S)*t)*(1-h_D*(1-h_S)*c_trait - h_S*c_sup)*x01*x10
        + (1-c_sup)*x01*x11
    )

    f_x10 = (
        (1+h_D*t)*(1-h_D*c_trait)*x00*x10
        + 0.5*(1+h_D*(1-h_S)*t)*(1-h_D*(1-h_S)*c_trait - h_S*c_sup)*x00*x11
        + 0.5*(1+h_D*(1-h_S)*t)*(1-h_D*(1-h_S)*c_trait - h_S*c_sup)*x01*x10
        + (1-c_trait)*x10**2
        + (1-(1-h_S)*c_trait - h_S*c_sup)*x10*x11
    )

    f_x11 = (
        0.5*(1+h_D*(1-h_S)*t)*(1-h_D*(1-h_S)*c_trait - h_S*c_sup)*x00*x11
        + 0.5*(1+h_D*(1-h_S)*t)*(1-h_D*(1-h_S)*c_trait - h_S*c_sup)*x01*x10
        + (1-c_sup)*x01*x11
        + (1-(1-h_S)*c_trait - h_S*c_sup)*x10*x11
        + (1-c_sup)*x11**2
    )

    w = f_x00 + f_x01 + f_x10 + f_x11
    return f_x00/w, f_x01/w, f_x10/w, f_x11/w


t = 0.87 * k
c_trait = 0.9 * k**1.5

x00, x01, x10, x11 = 1 - 2e-3, 0, 1e-3, 1e-3

x00_list, x01_list, x10_list, x11_list = [], [], [], []

for gen in range(max_generations):
    x00_list.append(x00)
    x01_list.append(x01)
    x10_list.append(x10)
    x11_list.append(x11)

    x00, x01, x10, x11 = next_gen(
        x00, x01, x10, x11,
        h_D, h_S, t, c_trait, c_sup
    )

plt.figure(figsize=(10, 5))
plt.plot(x00_list, label=r"$x_{00}$")
plt.plot(x01_list, label=r"$x_{01}$")
plt.plot(x10_list, label=r"$x_{10}$")
plt.plot(x11_list, label=r"$x_{11}$")

plt.title(fr"$h_D={h_D},\ h_S={h_S},\ k={k}$")
plt.xlabel("Generation")
plt.ylabel("Frequency")
plt.ylim(0, 1.05)
plt.legend()
plt.tight_layout()
plt.show()