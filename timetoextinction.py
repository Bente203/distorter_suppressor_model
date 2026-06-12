# -*- coding: utf-8 -*-
"""
Created on Sun Mar 22 22:13:23 2026

@author: Bente
"""

import numpy as np
import matplotlib.pyplot as plt


c_sup = 0.15
h_S = 0
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

#points until where distorter can invade
#marker_points = {
 #   0: 0.93,
 #   0.2: 0.73,
 #   0.4: 0.62, 
 #   0.6: 0.56,
 #   0.8: 0.51,
 #   1.0: 0.47
#}

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

   # k_dot = marker_points[h_D]
   # idx = np.argmin(np.abs(k_values - k_dot))

   # plt.scatter(k_values[idx], yvals[idx],
     #           color=line.get_color(), s=25, zorder=5)

plt.xlabel("Trait distorter strength (k)")
plt.ylabel("log10 time to extinction")


plt.axvline(x=0.3 , color='lightgreen', linestyle= '--', label= 'suppressor invasion')

plt.legend(fontsize = 5)
plt.show()


#%% 
#frequencies of different genotypes at values of k


import numpy as np
import matplotlib.pyplot as plt

c_sup = 0.15
h_S = 0.5

h_D=0.6

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
    
    return (
        f_x00/w,
        f_x01/w,
        f_x10/w,
        f_x11/w
    )

k_values = np.linspace(0,1,100)

plt.figure(figsize=(10,5))

distorter_final = []
suppressor_final = []

for k in k_values:

    t = 0.87 * k
    c_trait = 0.9 * k**1.5

    x00, x01, x10, x11 = 1 - 2e-6, 0, 1e-6, 1e-6 #0.97, 0.01, 0.01, 0.01

    for gen in range(max_generations):
        x00, x01, x10, x11 = next_gen(x00, x01, x10, x11)

        if x10 < threshold:
            break

    D = x10 
    S = x01

    distorter_final.append(D)
    suppressor_final.append(S)

plt.plot(k_values, distorter_final, label=f"D, hD={h_D}")
plt.plot(k_values, suppressor_final, linestyle="--", label=f"S, hS={h_S}")

plt.xlabel("k")
plt.ylabel("Final frequency")
plt.legend()
plt.title("Final frequencies of distorter and suppressor")
plt.show()
#%%

import numpy as np
import matplotlib.pyplot as plt

c_sup = 0.15
h_S = 0.5
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

#%%

import numpy as np

c_sup = 0.15
h_S = 0.5
h_D = 0.4

max_generations = 20000

def next_gen(x00, x01, x10, x11, t, c_trait):
    f_x00 = (
        x00**2 +
        x00*x01 + 
        (1-h_D*t)*(1-h_D*c_trait)*x00*x10 +
        1/2*(1-h_D*(1-h_S)*t)*(1-h_D*(1-h_S)*c_trait - h_S*c_sup)*x00*x11 +
        1/2*(1-h_D*(1-h_S)*t)*(1-h_D*(1-h_S)*c_trait - h_S*c_sup)*x01*x10 
    )
    
    f_x01 = (
        x00*x01 +
        1/2*(1-h_D*(1-h_S)*t)*(1-h_D*(1-h_S)*c_trait - h_S*c_sup)*x00*x11 +
        x01**2 +
        1/2*(1-h_D*(1-h_S)*t)*(1-h_D*(1-h_S)*c_trait - h_S*c_sup)*x01*x10 +
        (1-c_sup)*x01*x11
    )
    
    f_x10 = (
        (1+h_D*t)*(1-h_D*c_trait)*x00*x10 +
        1/2*(1+h_D*(1-h_S)*t)*(1-h_D*(1-h_S)*c_trait - h_S*c_sup)*x00*x11 +
        1/2*(1+h_D*(1-h_S)*t)*(1-h_D*(1-h_S)*c_trait - h_S*c_sup)*x01*x10 + 
        (1-c_trait)*x10**2 +
        (1-(1-h_S)*c_trait - h_S*c_sup)*x10*x11   
    )
    
    f_x11 = (
        1/2*(1+h_D*(1-h_S)*t)*(1-h_D*(1-h_S)*c_trait - h_S*c_sup)*x00*x11 +
        1/2*(1+h_D*(1-h_S)*t)*(1-h_D*(1-h_S)*c_trait - h_S*c_sup)*x01*x10 +
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


k = 0.35

t = 0.87 * k
c_trait = 0.9 * k**1.5

x00, x01, x10, x11 = 1 - 2e-3, 0, 1e-3, 1e-3

for gen in range(max_generations):
    x00, x01, x10, x11 = next_gen(x00, x01, x10, x11, t, c_trait)

print(f"k = {k}")
print(f"x00 = {x00:.10f}")
print(f"x01 = {x01:.10f}")
print(f"x10 = {x10:.10f}")
print(f"x11 = {x11:.10f}")
print(f"som = {x00 + x01 + x10 + x11:.10f}")