#!/usr/bin/env python3
import numpy as np
import matplotlib.pyplot as plt
import time
from random import choice

from numpy.linalg import inv
from scipy.optimize import linprog
from scipy import interpolate
import scipy.integrate as integrate


def plot_profil_tab(pt_mes, alt, tab):

    plt.plot(pt_mes, alt, 'k', label="Profil du terrain")
    for _ in range(100):
      profil = choice(tab)
      plt.plot(pt_mes,profil,'grey')
    plt.plot(pt_mes,profil,'grey',label="Profil du terrain, itération d'Uzawa")

    xnew, profil_liss = lisse_moi(pt_mes, profil)

    plt.plot(xnew, profil_liss, 'b', label="Route optimale lissée")

    plt.xlabel("Points de mesures (en metres)")
    plt.ylabel("Altitudes (en metres)")
    plt.title("Tracés du profil du terrain et de la route optimale associée")
    plt.legend()
    plt.show()


def plot_profil(pt_mes, alt, profil):

    plt.plot(pt_mes, alt, 'k', label="Profil du terrain")
    plt.plot(pt_mes, profil, 'c', label="Route optimale")

    xnew, profil_liss = lisse_moi(pt_mes, profil)

    plt.plot(xnew, profil_liss, 'b', label="Route optimale lissée")

    plt.xlabel("Points de mesures (en metres)")
    plt.ylabel("Altitudes (en metres)")
    plt.title("Tracés du profil du terrain et de la route optimale associée")
    plt.legend()
    plt.show()

def pente_max(pt_mes, alt):

    lenth = len(pt_mes)

    penteMax = abs(alt[0]-alt[1])/abs(pt_mes[0]-pt_mes[1])

    for i in range(1,lenth):
        penteTemp = abs(alt[i-1]-alt[i])/abs(pt_mes[i-1]-pt_mes[i])
        if penteMax<penteTemp:
            penteMax=penteTemp

    return penteMax

  
def lisse_moi(pt_mes, profil, xcoeff=20):
    """lisse le profil optimal en utilisant Spline
    param:
    :pt_mes: liste des points
    :profil: profil optimal à lisser
    :xcoeff: nouveau nombre de points xcoeff*len(pt_mes)
    return: array 
    """
    tck = interpolate.splrep(pt_mes, profil, s=len(pt_mes)-np.sqrt(2*len(pt_mes)))
    new_nb_pt = xcoeff*len(pt_mes)
    xnew = np.linspace(pt_mes[0], pt_mes[-1], new_nb_pt)
    profil_liss = interpolate.splev(xnew, tck)

    return xnew, profil_liss

def lire_profil(fichier):

    with open(fichier, 'r') as data:
        pt_mes = []
        alt = []
        for line in data:
            p = line.split()
            pt_mes.append(float(p[0]))
            alt.append(float(p[1]))
    return pt_mes, alt

lon,alt=lire_profil("profil.txt")

def def_matrice_quad(pt_mes, alt, h, penteMax):

    ATemp = [h/4]+(len(pt_mes)-2)*[h/2]+[h/4]
    A = np.diag(ATemp)

    inv_A = inv(A)

    b = np.array((2*(len(pt_mes)-1))*[penteMax*h])

    C = np.zeros((2*(len(pt_mes)-1), len(pt_mes)))
    for i in range(len(pt_mes)-1):
        C[i][i] = -1
        C[i][i+1] = 1
        C[i+(len(pt_mes)-1)][i] = 1
        C[i+(len(pt_mes)-1)][i+1] = -1

    return A, inv_A, b, C


penteMax=0.1

h=1
A,inv_A,b,C=def_matrice_quad(lon,alt,h,penteMax)


    
def uzawa(A,inv_A,alt,C,b,rho=0.2,epsi=0.01,maxi=20000) :
   n=len(alt)
   l=len(b)
   tab = []
   #print(l)

   alt = np.array(alt)
   
   xold = np.ones(n).T
   zold = np.ones(l).T
   
   x = alt-(1/2)*inv_A.dot(((C.T).dot(zold)))

   print(n,len(x))


   b= np.array(b)
   print(len(zold+rho*(C.dot(x.T)-b)))

   z = np.maximum(np.zeros((l)),zold+rho*(C.dot(x.T)-b))

   k = 0
   while (k <= maxi and np.linalg.norm(xold - x) > epsi) :
       tab.append(x)
       k += 1
       xold = x
       x=alt-(1/2)*inv_A.dot((C.T).dot(z))

       z=np.maximum(np.zeros(l),z+rho*(C.dot(x.T)-b))

   if k>= maxi:
       print("borne atteinte")
   else:
       print(k)
       
   
   return x,tab,k
   
def evolv_alpha(fichier="", methode='linprog', alpha=0.1, pas=5):
    """affiche le trace de la route optimale en fonction de la pente maximale
    param:
    :fichier:route a optimiser (point mesure/altitude associee - separees par des espaces)
    :methode: linprog (utilise scipy.linprog) ou maison (utilise un simplexe maison)
    :pas: affiche le trace de la route optimale trouvee pour :pas: pentes maximales
    """

    pt_mes, alt = lire_profil('profil.txt')
    plt.plot(pt_mes, alt, 'k')

    L = pt_mes[-1]
    n = len(pt_mes)
    h = L/(n-1)

    penteLoc = pente_max(pt_mes, alt) #pente maximale de la route
    alpha = 0.1 #alpha = 10% fixe
    evantail_pente = np.linspace(penteLoc-(penteLoc+alpha)/2, alpha, pas)
    pente_testee = []

    change_color=1/(pas+1) #intensifie le bleu des tracés au cours des iterations

    for pente in evantail_pente:

        A,inv_A,b,C=def_matrice_quad(lon,alt,h,pente+alpha) #matrices

        bnds = [(None, None)] * n + [(0, None)] * n #bornes

        profil=uzawa(A,inv_A,alt,C,b)

        pente_testee.append('alpha='+str(pente)) #legend

        plt.plot(pt_mes, profil, color=(0,0,change_color)) #affichage route optimale trouvee
        change_color+=1/(pas+1)

    #affichage route optimale lissée
    xnew, profil_liss = lisse_moi(pt_mes, profil)
    plt.plot(xnew, profil_liss, 'r')


    leg=np.concatenate((["Profil du terrain"], pente_testee, ["Route optimale lissée"]))
    plt.legend(leg)
    plt.xlabel("Points de mesures (en metres)")
    plt.ylabel("Altitudes (en metres)")
    plt.title("Evolution du tracé de la route optimale en fonction de la pente maximale alpha")
    plt.show()


#evolv_alpha()
#b,tab,k=uzawa(A,inv_A,alt,C,b)


def plot_alpha():
  alphas = np.linspace(0,0.10,num=100)
  tab = []
  for alpha in alphas:

    h=1
    A,inv_A,b,C=def_matrice_quad(lon,alt,h,alpha)
    _,_,k = uzawa(A,inv_A,alt,C,b)
    tab.append(k)

  plt.plot(alphas, tab)

  plt.xlabel("Alpha")
  plt.ylabel("Nombre d'iterations $k$")
  plt.title("Nombre d'itérations de l'algorithme d'Uzawa en fonction de $\\alpha$")
  plt.legend()
  plt.show()






#plot_profil_tab(lon,alt,tab)