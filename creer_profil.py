# -*- coding: utf-8 -*-
"""
Created on Sun Oct  7 16:50:59 2018

@author: Xavier
"""
import numpy as np
import random
def generate_profil():
    dist=np.linspace(200,5000,100)
    a=random.random()*200
    delta=1.2
    alt=[a]
    for i in range(1000-1):
        alt.append(alt[-1]+random.random()*delta*(-1)**random.randint(1,2))
        
    print(dist)
    with open("profil3.txt",'w') as file:
        for i in range(100):
            file.write("    {}    {}\n".format(dist[i],alt[i]))
    
generate_profil()