import numpy as np
import matplotlib.pyplot as plt
import time

from numpy.linalg import inv
from scipy.optimize import linprog
from scipy import interpolate
import scipy.integrate as integrate

from simplexeX import simplex
from jacobo import uzawa

def lire_profil(fichier):

    with open(fichier, 'r') as data:
        pt_mes = []
        alt = []
        for line in data:
            p = line.split()
            pt_mes.append(float(p[0]))
            alt.append(float(p[1]))
    return pt_mes, alt

def plot_profil(pt_mes, alt, profil):

    plt.plot(pt_mes, alt, 'k', label="Profil du terrain")
    plt.plot(pt_mes, profil, 'c', label="Route optimale")

    xnew, profil_liss = lisse_moi(pt_mes, profil)

    plt.plot(xnew, profil_liss, 'b', label="Route optimale lissée")

    plt.xlabel("Points de mesures (en metres)")
    plt.ylabel("Altitudes (en metres)")
    plt.title("Tracés du profil du terrain et de la route optimale associée")
    plt.legend()

def plot_couts(fichier, profil):

    pt_mes, alt = lire_profil(fichier)

    total_terre=[alt[i] - profil[i] for i in range(len(alt))]

    plt.plot(pt_mes, total_terre)
    plt.xlabel("Points de mesures (en metres)")
    plt.ylabel("Quantité de terre déplacée (en metres)")
    plt.title("Quantité de terre déplacée le long du terrain")
    plt.legend()
    plt.show()

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

def pente_max(pt_mes, alt):

    lenth = len(pt_mes)

    penteMax = abs(alt[0]-alt[1])/abs(pt_mes[0]-pt_mes[1])

    for i in range(1,lenth):
        penteTemp = abs(alt[i-1]-alt[i])/abs(pt_mes[i-1]-pt_mes[i])
        if penteMax<penteTemp:
            penteMax=penteTemp

    return penteMax


    #CAS LINEAIRE

def def_matrice(pt_mes, alt, h, penteMax):
    """Soit le probleme d'optimisation min f^t * v
    s.c C * v - d <= 0
    return: f, C, d
    """
    
    n = len(pt_mes)

    f = np.concatenate((np.zeros(n), [1/2], np.ones(n-2), [1/2]))
    f = np.multiply(h, f)

    C = np.zeros( (2*(2*n-1), 2*n) )
    for i in range(n-1):
        C[i][i] = -1
        C[i][i+1] = 1
        C[i+(n-1)][i] = 1
        C[i+(n-1)][i+1] = -1

        C[i+2*(n-1)][i] = 1
        C[i+3*n-2][i] = -1

        C[i+2*(n-1)][i+n] = -1/h
        C[i+3*n-2][i+n] = -1/h

    C[3*(n-1)][n-1] = 1
    C[(n-1)+3*n-2][n-1] = -1
    C[3*(n-1)][2*n-1] = -1/h
    C[(n-1)+3*n-2][2*n-1] = -1/h

    d = np.concatenate((np.array((2*(n-1))*[penteMax*h]), alt, np.multiply(alt,-1)))    

    return f, C, d

def find_profil_lin(fichier, methode='linprog', alpha = 0.1,verbose=False):
    """affiche le trace de la route et le profil optimal associe
    param:
    :fichier: route a optimiser (point mesure/altitude associee - separees par des espaces)
    :methode: linprog (utilise scipy.linprog) ou maison (utilise un simplexe maison)
    :verbose: affiche la progression du programme ainsi que les matrices generees
    """
    time_sleep = 0.5

    if verbose:
        print('\n')
        print("Lecture du fichier...\n")
        time.sleep(time_sleep)
    pt_mes, alt = lire_profil(fichier)

    L = pt_mes[-1]
    n = len(pt_mes)
    h = L/(n-1)

    penteLoc = pente_max(pt_mes, alt) #pente maximale du terrain
    if verbose:
        print('\n')
        print("Calcul de la pente locale maximale...\n")
        print ("Pente locale maximale trouvée : ", penteLoc)
        time.sleep(time_sleep)

    alpha = 0.1 #alpha = 10% fixe

    f, C, d = def_matrice(pt_mes, alt, h, alpha) #matrices
    if verbose:
        print('\n')
        print ("Le probleme d'optimisation considere est :\n")
        print ("min f^t * v \ns.c C * v - d <= 0 \n")
        time.sleep(time_sleep)

    if verbose:
        print('\n')
        print ("Calcul du profil optimal associe...\n")
        time.sleep(time_sleep)
    if methode=='linprog':
        bnds = [(None, None)] * n + [(None, None)] * n #bornes
        #bnds=tuple(bnds)
        #print(bnds)
        res = linprog(f, A_ub=C,b_ub=d,bounds=bnds,options={'maxiter':1700,'tol':10e-9})
        print(res["nit"])
        profil = res['x'][:n] #profil retournee par scipy.linprog
    if methode=='maison':
        res = simplex(f, C, d, verbose)
        print(len(res))
        profil = res[:n] #profil retournee par le simplexe maison
    print("min={}".format(sum([abs(alt[k]-profil[k]) for k in range(len(alt))])))
    if verbose:
        print('\n')
        print ("Route", alt)
        print("Profil associé", profil)

    plot_profil(pt_mes, alt, profil) #affichage

    plt.show()

    return profil

def evolv_alpha(fichier, methode='linprog', alpha=0.1, pas=5):
    """affiche le trace de la route optimale en fonction de la pente maximale
    param:
    :fichier:route a optimiser (point mesure/altitude associee - separees par des espaces)
    :methode: linprog (utilise scipy.linprog) ou maison (utilise un simplexe maison)
    :pas: affiche le trace de la route optimale trouvee pour :pas: pentes maximales
    """

    pt_mes, alt = lire_profil(fichier)
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

        f, C, d = def_matrice(pt_mes, alt, h, pente) #matrices

        bnds = [(None, None)] * n + [(0, None)] * n #bornes

        if methode=='linprog':
            res = linprog(f, A_ub=C,b_ub=d,bounds=bnds,options={'maxiter':10000,'tol':10e-9})
            print(res)
            profil = res['x'][:n] #profil retournee par scipy.linprog
        if methode=='maison':
            res = simplex(f, C, d)
            profil = res[:n] #profil retournee par le simplexe maison

        pente_testee.append('alpha='+str(round(pente,2))) #legend

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


    #CAS QUADRATIQUE

def def_matrice_quad(pt_mes, alt, h, penteMax):
    """Soit le probleme d'optimisation min <A(U-G), (U-G)>
    s.c C * U - d <= 0
    return: A, inv_A, C, d
    """

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

def main():
    
    #1fichier = input("Veuillez entrer le nom du fichier :\n> ")
    fichier='profil3.txt'

    #choix = int(input("Cas linéaire (0) ou cas quadratique (1) :\n> "))
    choix=0

    if not choix:   
    #CAS LINEAIRE

        #methode = int(input("Utiliser scipy.linprog (0) ou simplexe maison (1) :\n> "))
        methode=1
        if not methode:
            profil = find_profil_lin(fichier, verbose=False) #scipy.linprog

            """plot = input("Quantité de terre déplacée le long du terrain? (o/n)\n> ")
            if plot=='o':
                plot_couts(fichier, profil)

            plot = input("Suivre l'evolution du tracés en fonction de la pente maximale? (o/n)\n> ")
            if plot=='o':
                evolv_alpha(fichier)"""

        else:
            print("==Veuillez patienter, le programme est un peu long... (environ 1000 itérations)==\n")
            time.sleep(0.5)
            profil = find_profil_lin(fichier, methode='maison', verbose=True) #simplexe maison

            plot = input("Quantité de terre déplacée le long du terrain? (o/n)\n> ")
            if plot=='o':
                plot_couts(fichier, profil)

            plot_evol = input("Suivre l'evolution du tracés en fonction de la pente maximale? (o/n)\n==RISQUE DE PRENDRE BEAUCOUP DE TEMPS==\n> ")
            if plot_evol=='o':
                evolv_alpha(fichier, methode='maison')

    else:
    #CAS QUADRATIQUE

        pt_mes, alt = lire_profil(fichier)

        L = pt_mes[-1]
        n = len(pt_mes)
        h = L/(n-1)
        print("Cas partiellement traité...\n")
        print("Calcul des matrices...\n")
        A, inv_A, b, C = def_matrice_quad(pt_mes, alt, h, penteMax=0.1)
        x=uzawa(A,inv_A,alt,C,b)
        btemp=x[0][0]
        #print(btemp)
        #print(len(btemp))
        #print(len(pt_mes))
        #print(alt)
        plot_profil(pt_mes,alt,btemp.T)
if __name__=="__main__":
    main()