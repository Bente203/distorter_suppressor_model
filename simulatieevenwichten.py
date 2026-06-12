# -*- coding: utf-8 -*-
"""
Created on Sun Mar 15 13:50:52 2026

@author: Bbent
"""
import numpy as np
import matplotlib.pyplot as plt


c_sup = 0.15
h_Dvalues = [0, 0.2, 0.3, 0.4, 0.6, 0.8, 1.0]
h_S = 0

iterations = 20000

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
for h_D in h_Dvalues:
    distortion_eq = []
    
        
    
    for k in k_values:
    
        t = 0.87 * k
        c_trait = 0.9 * k**1.5
    
        x00, x01, x10, x11 = 1 - 2e-3, 0, 1e-3, 1e-3
        
        for i in range(iterations):
            x00, x01, x10, x11 = next_gen(x00, x01, x10, x11)
        
    
        distortion_eq.append(x10+x11)

    plt.plot(k_values, distortion_eq, label=f"h_D={h_D}")
max_index = np.argmax(distortion_eq)
k_max = k_values[max_index]

plt.xlabel("Trait distorter strength (k)")
plt.ylabel("Distorter frequency (x10+x11)") #(x10*k + x11*k)
plt.legend()
plt.show()


print(x00, x01, x10, x11)
print(k_max)

    

#%%

import numpy as np
import matplotlib.pyplot as plt

c_sup = 0.15
h_S = 1
h_D = 1

max_generations = 20000

def next_gen(x00, x01, x10, x11, t, c_trait):
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
    
    return (
        f_x00/w,
        f_x01/w,
        f_x10/w,
        f_x11/w
    )

k_values = np.linspace(0.001, 1, 400)

x00_final = []
x01_final = []
x10_final = []
x11_final = []

for k in k_values:
    t = 0.87 * k
    c_trait = 0.9 * k**1.5

    x00, x01, x10, x11 = 1 - 2e-3, 0, 1e-3, 1e-3
    
    for gen in range(max_generations):
        x00, x01, x10, x11 = next_gen(x00, x01, x10, x11, t, c_trait)

    x00_final.append(x00)
    x01_final.append(x01)
    x10_final.append(x10)
    x11_final.append(x11)

fig, ax = plt.subplots(2, 2, figsize=(11, 7), sharex=True)

ax[0,0].plot(k_values, x00_final, color='black')
ax[0,0].set_title(r"$x_{00}$")

ax[0,1].plot(k_values, x01_final, color='orange')
ax[0,1].set_title(r"$x_{01}$")

ax[1,0].plot(k_values, x10_final, color='blue')
ax[1,0].set_title(r"$x_{10}$")

ax[1,1].plot(k_values, x11_final, color='red')
ax[1,1].set_title(r"$x_{11}$")

for a in ax.ravel():
    a.set_ylim(-0.02, 1.02)
    a.grid(alpha=0.3)

ax[1,0].set_xlabel("k")
ax[1,1].set_xlabel("k")
ax[0,0].set_ylabel("Final frequency")
ax[1,0].set_ylabel("Final frequency")

plt.suptitle(fr"All final gamete frequencies ($h_D={h_D}$, $h_S={h_S}$)")
plt.tight_layout()
plt.show()












    