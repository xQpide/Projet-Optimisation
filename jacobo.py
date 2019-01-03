#!/usr/bin/env python3
import numpy as np
import matplotlib.pyplot as plt
import time

from numpy.linalg import inv
from scipy.optimize import linprog
from scipy import interpolate
import scipy.integrate as integrate



    
def uzawa(A,inv_A,alt,C,b,rho=0.02,epsi=0.001,maxi=50000) :
    
   n=len(alt)
   l=len(b)
   alt=np.array(alt)
   xold = np.matrix(np.ones(n)).T
   zold = np.matrix(np.ones(l)).T
   print(xold)
   x =alt-(1/2)*np.dot(inv_A,np.dot(C.T,zold))
   print(x)
   b= np.matrix(b).T
   #print(C,'\n',b)
   #print(zold+rho*(C.dot(x)-b))
   #print(np.zeros((l,l)).T)
   z = np.maximum(np.zeros((185,l)).T,zold+rho *(C.dot(x)-b))
   print(len(z))
   k = 0
   while (k <= maxi and np.linalg.norm(xold - x) > epsi) :
       #♥print("allo{}".format(k))
       print(alt-(1/2)*np.dot(inv_A,np.dot(C.T,z)))
       k =k+ 1
       xold = x
       x=alt-(1/2)*np.dot(inv_A,np.dot(C.T,z))
       z=np.maximum(np.zeros((185,l)).T,z+rho*(C.dot(x)-b))
       print(np.linalg.norm(xold - x))
       
   
   return x
   

if __name__=="__main__":
    b=uzawa(A,inv_A,alt,C,b)
    
    plot_profil(lon,alt,b)