# -*- coding: utf-8 -*-
"""
Created on Sun Mar 15 13:50:52 2026

@author: Bbent
"""
import numpy as np
import matplotlib.pyplot as plt


c_sup = 0.15
h_Dvalues = [0, 0.2, 0.3, 0.4, 0.6, 0.8, 1.0]
h_S = 0.5

iterations = 100000


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
    
        x00, x01, x10, x11 = 1 - 2e-6, 0, 1e-6, 1e-6
        
        for i in range(iterations):
            x00, x01, x10, x11 = next_gen(x00, x01, x10, x11)
        
    
        # equilibrium distortion opslaan
        distortion_eq.append(x10+x11) #(x10*k + x11*k)

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

h_D = 0.2
c_sup = 0.15

k_values = np.linspace(0.001, 1, 1000)

p3_values =[]
distortion_p3 =[]

for k in k_values:
    t = 0.87 * k
    c_trait = 0.9 * k**1.5

    p3 = h_D * (c_trait * (1 + h_D * t) - t) / (c_trait * (2*h_D - 1))

    p3_values.append(p3)
    distortion_p3.append(p3)

p3_values = np.array(p3_values)
distortion_p3 = np.array(distortion_p3)

mask = (p3_values >= 0) & (p3_values <= 1)


#plt.figure(figsize=(8,5))
plt.plot(k_values[mask], distortion_p3[mask], color='darkblue', label="$p_3^*$")

plt.xlabel("Trait distorter strength (k)")
plt.ylabel("Predicted distortion at internal equilibrium $p_3^*$")
plt.legend()
plt.tight_layout()
plt.show()

    
#%%

import numpy as np
import matplotlib.pyplot as plt

# parameters
c_sup = 0.15
h_S = 1
h_D = 0.4       
k = 0.5           
t = 0.87 * k
c_trait = 0.9 * k**1.5

iterations = 3000

def next_gen(x00, x01, x10, x11, h_D, h_S, t, c_trait, c_sup):
    f_x00 = (
        x00**2 +
        x00*x01 +
        (1-h_D*t)*(1-h_D*c_trait)*x00*x10 +
        0.5*(1-h_D*(1-h_S)*t)*(1-h_D*(1-h_S)*c_trait - h_S*c_sup)*x00*x11 +
        0.5*(1-h_D*(1-h_S)*t)*(1-h_D*(1-h_S)*c_trait - h_S*c_sup)*x01*x10
    )

    f_x01 = (
        x00*x01 +
        0.5*(1-h_D*(1-h_S)*t)*(1-h_D*(1-h_S)*c_trait - h_S*c_sup)*x00*x11 +
        x01**2 +
        0.5*(1-h_D*(1-h_S)*t)*(1-h_D*(1-h_S)*c_trait - h_S*c_sup)*x01*x10 +
        (1-c_sup)*x01*x11
    )

    f_x10 = (
        (1+h_D*t)*(1-h_D*c_trait)*x00*x10 +
        0.5*(1+h_D*(1-h_S)*t)*(1-h_D*(1-h_S)*c_trait - h_S*c_sup)*x00*x11 +
        0.5*(1+h_D*(1-h_S)*t)*(1-h_D*(1-h_S)*c_trait - h_S*c_sup)*x01*x10 +
        (1-c_trait)*x10**2 +
        (1-(1-h_S)*c_trait - h_S*c_sup)*x10*x11
    )

    f_x11 = (
        0.5*(1+h_D*(1-h_S)*t)*(1-h_D*(1-h_S)*c_trait - h_S*c_sup)*x00*x11 +
        0.5*(1+h_D*(1-h_S)*t)*(1-h_D*(1-h_S)*c_trait - h_S*c_sup)*x01*x10 +
        (1-c_sup)*x01*x11 +
        (1-(1-h_S)*c_trait - h_S*c_sup)*x10*x11 +
        (1-c_sup)*x11**2
    )

    w = f_x00 + f_x01 + f_x10 + f_x11

    return f_x00 / w, f_x01 / w, f_x10 / w, f_x11 / w


# beginwaarden
x00, x01, x10, x11 = 1 - 2e-6, 0, 1e-6, 1e-6

# opslag
time = []
x00_values = []
x01_values = []
x10_values = []
x11_values = []

for i in range(iterations):
    time.append(i)
    x00_values.append(x00)
    x01_values.append(x01)
    x10_values.append(x10)
    x11_values.append(x11)

    x00, x01, x10, x11 = next_gen(x00, x01, x10, x11, h_D, h_S, t, c_trait, c_sup)

# plot
plt.figure(figsize=(10, 6))
plt.plot(time, x00_values, label=r"$x_{00}$")
plt.plot(time, x01_values, label=r"$x_{01}$")
plt.plot(time, x10_values, label=r"$x_{10}$")
plt.plot(time, x11_values, label=r"$x_{11}$")

plt.xlabel("Generations")
plt.ylabel("Frequency")
plt.title(f"Time dynamics for h_D={h_D}, h_S={h_S}, k={k}")
plt.legend()
plt.tight_layout()
plt.show()


















    