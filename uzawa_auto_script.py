import numpy as np
from scipy.optimize import minimize

# données
n = 2
m = 1
C = np.array([1,1])  # un matrice m*n
b = np.array([-1])    # un vecteur de R^m

def l_exemple(u,a): # le Lagrengien, il prend deux parametre u E R^n et a E R^m
	# a changer selon le probleme
	return np.linalg.norm(u)**2 + a * (u[0] + u[1] - 1)

def cl(x):  # a changer selon le probleme
	return np.dot(C, x) - b


def mini_python(l,x0, z):  # pour calculer le min dans etape (i)
	def l_zk(x_):
		return l(x_, z)
	res = minimize(l_zk, x0, method='nelder-mead')
	return res.x

def uzawa(l, p, e, kmax):
   """
    l : le Lagrangien, fonction qui prend deux parametre x un vecteur de dimension n z de dimension m
    p : pas
    e : tolerence d'erreur
    kmax : nb_max d'teration
   """
   x = np.zeros(n)
   z = np.ones(m)
   xk = mini_python(l,x, z)
   zk = max(0, z + p * cl(xk))
   k = 1
	
   while (k <= kmax and np.linalg.norm(x - xk) > e):
       x = xk
       xk = mini_python(l,xk, zk)
       zk = max(0, zk + p * cl(xk))
       k += 1
       
   return xk, zk


ress = uzawa(l_exemple, 1, 0.001, 22)
print(ress)
