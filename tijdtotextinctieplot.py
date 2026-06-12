# -*- coding: utf-8 -*-
"""
Created on Sun Mar 22 22:13:23 2026

@author: Bente
"""

import numpy as np
import matplotlib.pyplot as plt


c_sup = 0.15
h_S = 0.5
hD_values = [0, 0.2, 0.3,0.4, 0.6, 0.8, 1.0]

threshold = 1e-6
max_generations = 20000
def next_gen(x00, x01, x10, x11):
    
    f_x00 = (
        x00**2 +
        x00*x01 + 
        (1-h_D * t)*(1-h_D*c_trait)*x00*x10 +
        1/2 * (1-h_D*(1-h_S)*t)*(1-h_D*(1-h_S)*c_trait - h_S*c_sup)*x00*x11 +
        1/2 * (1-h_D*(1-h_S)*t)*(1-h_D*(1-h_S)*c_trait - h_S*c_sup)*x01*x10 
        )
    
    f_x01 = (
        x00*x01 +
        1/2 * (1-h_D*(1-h_S)*t)*(1-h_D*(1-h_S)*c_trait - h_S*c_sup)*x00*x11 +
        x01**2 +
        1/2 * (1-h_D*(1-h_S)*t)*(1-h_D*(1-h_S)*c_trait - h_S*c_sup)*x01*x10 +
        (1-c_sup)*x01*x11
        )
    
    f_x10 = (
        (1+h_D * t)*(1-h_D*c_trait)*x00*x10 +
        1/2 * (1+h_D*(1-h_S)*t)*(1-h_D*(1-h_S)*c_trait - h_S*c_sup)*x00*x11 +
        1/2 * (1+h_D*(1-h_S)*t)*(1-h_D*(1-h_S)*c_trait - h_S*c_sup)*x01*x10 + 
        (1-c_trait)*x10**2 +
        (1-(1-h_S)*c_trait - h_S*c_sup)*x10*x11   
        )
    
    f_x11 = (
        1/2 * (1+h_D*(1-h_S)*t)*(1-h_D*(1-h_S)*c_trait - h_S*c_sup)*x00*x11 +
        1/2 * (1+h_D*(1-h_S)*t)*(1-h_D*(1-h_S)*c_trait - h_S*c_sup)*x01*x10 +
        (1-c_sup)*x01*x11 +
        (1-(1-h_S)*c_trait - h_S*c_sup)*x10*x11 +
        (1-c_sup)*x11**2
        ) 
    w = f_x00 + f_x01 + f_x10 + f_x11
    
    x00_new = f_x00 / w
    x01_new = f_x01 / w
    x10_new = f_x10 / w
    x11_new = f_x11 / w

    return (
        x00_new,
        x01_new,
        x10_new,
        x11_new       
        
        )


k_values = np.linspace(0,1,100)

for h_D in hD_values:

    time_to_extinction = []

    for k in k_values:

        t = 0.87 * k
        c_trait = 0.9 * k**1.5


        x00, x01, x10, x11 = 1 - 2e-3, 0, 1e-3, 1e-3

        for gen in range(max_generations):
            x00, x01, x10, x11 = next_gen(x00, x01, x10, x11)

            if x10 + x11 < threshold:
                break

        time_to_extinction.append(gen)

    yvals = np.log10(time_to_extinction)
    line, = plt.plot(k_values, yvals, label=f"h_D={h_D}")

plt.xlabel("Trait distorter strength (k)")
plt.ylabel("log10 time to extinction")
plt.legend(fontsize = 5)
plt.show()

