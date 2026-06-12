# -*- coding: utf-8 -*-
"""
Created on Sun Mar 15 13:50:52 2026

@author: Bbent
"""
import numpy as np
import matplotlib.pyplot as plt


c_sup = 0.15
h_Dvalues = [0, 0.2, 0.5, 0.8, 1.0]
h_S = 0

iterations = 10000


def next_gen(x00, x01, x10, x11):
    w = (
    x00**2
    + x01**2
    + (1 - c_trait) * x10**2
    + (1 - c_sup) * x11**2
    + 2 * x00 * x01
    + 2 * (1 - h_D * c_trait) * x00 * x10
    + 2 * (1 - h_D * (1 - h_S) * c_trait - h_S * c_sup) * x00 * x11
    + 2 * (1 - h_D * (1 - h_S) * c_trait - h_S * c_sup) * x10 * x01
    + 2 * (1 - c_sup) * x01 * x11
    + 2 * (1 + (1 - h_S) * t) * (1 - (1 - h_S) * c_trait - h_S * c_sup) * x10 * x11
)
    
    x00_new = (
        x00**2 +
        x00*x01 + 
        (1-h_D * t)*(1-h_D*c_trait)*x00*x10 +
        1/2 * (1-h_D*(1-h_S)*t)*(1-h_D*(1-h_S)*c_trait - h_S*c_sup)*x00*x11 +
        1/2 * (1-h_D*(1-h_S)*t)*(1-h_D*(1-h_S)*c_trait - h_S*c_sup)*x01*x10 
        ) /w
    
    x01_new = (
        x00*x01 +
        1/2 * (1-h_D*(1-h_S)*t)*(1-h_D*(1-h_S)*c_trait - h_S*c_sup)*x00*x11 +
        x01**2 +
        1/2 * (1-h_D*(1-h_S)*t)*(1-h_D*(1-h_S)*c_trait - h_S*c_sup)*x01*x10 +
        (1-c_sup)*x01*x11
        )/ w
    
    x10_new = (
        (1+h_D * t)*(1-h_D*c_trait)*x00*x10 +
        1/2 * (1+h_D*(1-h_S)*t)*(1-h_D*(1-h_S)*c_trait - h_S*c_sup)*x00*x11 +
        1/2 * (1+h_D*(1-h_S)*t)*(1-h_D*(1-h_S)*c_trait - h_S*c_sup)*x01*x10 + 
        (1-c_trait)*x10**2 +
        (1+(1-h_S)*t)*(1-(1-h_S)*c_trait - h_S*c_sup)*x10*x11   
        )/w
    
    x11_new = (
        1/2 * (1+h_D*(1-h_S)*t)*(1-h_D*(1-h_S)*c_trait - h_S*c_sup)*x00*x11 +
        1/2 * (1+h_D*(1-h_S)*t)*(1-h_D*(1-h_S)*c_trait - h_S*c_sup)*x01*x10 +
        (1-c_sup)*x01*x11 +
        (1+(1-h_S)*t)*(1-(1-h_S)*c_trait - h_S*c_sup)*x10*x11 +
        (1-c_sup)*x11**2
        )/w 
    
    return (
        x00_new,
        x01_new,
        x10_new,
        x11_new       
        
        )
k_values = np.linspace(0,1,100)
for h_D in h_Dvalues:
    distortion_eq = []
    
        
    
    for k in k_values:
    
        t = 0.87 * k
        c_trait = 0.9 * k**1.5
    
        x00, x01, x10, x11 = 1 - 2e-6, 0, 1e-6, 1e-6
        
        for i in range(iterations):
            x00, x01, x10, x11 = next_gen(x00, x01, x10, x11)
        
    
        # equilibrium distortion opslaan
        distortion_eq.append(x10*k + x11*k)

    plt.plot(k_values, distortion_eq, label=f"h_D={h_D}")
max_index = np.argmax(distortion_eq)
k_max = k_values[max_index]

plt.xlabel("Trait distorter strength (k)")
plt.ylabel("Individual trait distortion (x10*k + x11*k)")
plt.legend()
plt.show()

print(x00, x01, x10, x11)
print(k_max)

    
    
    