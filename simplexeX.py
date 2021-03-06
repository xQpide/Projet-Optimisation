import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import linprog
import time

def find_max(coeff):
    """ determine le variable entrante"""
    m=coeff[0]
    kstock=0
    for i in range(1,len(coeff)):
        if coeff[i]<m:
            m=coeff[i]
            kstock=i
    if m>=0:
        return None
    return kstock

def find_var_sortante(A,b,index):
    """ determine la variable sortante"""
    istart=0
    while istart<len(b)and A[istart][index]<=0:
        istart+=1
    if istart==len(b):
        return None
    m=b[istart]/A[istart][index]
    index_sortant=istart
    for i in range(istart,len(b)):
        if A[i][index]>0:
            if b[i]/A[i][index]<m:
                m=b[i]/A[i][index]
                index_sortant=i
    return index_sortant

def pivot(val_entrant,val_sortant,A,var_base,base,VE,coeff,b,objec):
    zero=np.zeros(len(A))
    """transformation du tableau par le pivot determiné par les variables entrantes et sortantes"""
    val_entrant=base.index(val_entrant)
    val_sortant=var_base.index(val_sortant)
    val_pivot=A[val_sortant][val_entrant]
    tmp1=A.T[val_entrant]
    tmp2=A[val_sortant]
    out=np.outer(tmp1,tmp2)

    Atemp=val_pivot*A-out
    
    zero[val_sortant]=Atemp[val_sortant][val_entrant]
    Atemp.T[val_entrant]=zero

    Atemp.T[val_entrant]=np.zeros(len(A))
    Atemp[val_sortant]=A[val_sortant]

    btemp=val_pivot*b-b[val_sortant]*tmp1
    btemp[val_sortant]=b[val_sortant]

    coefftemp=val_pivot*coeff-coeff[val_entrant]*tmp2 
    coefftemp[val_entrant]=0

    valtemp=b[val_sortant]/val_pivot
    objectemp=objec-coeff[val_entrant]*valtemp
    return Atemp/(val_pivot),coefftemp/(val_pivot),btemp/(val_pivot),objectemp

def init_simplex(f,A,b):
    cpt=0
    list_var=[]
    for i in range(len(f)):
        list_var.append("X{}".format(i))
    var_exces=[]
    var_artificielle=[]
    index_const_neg=[]
    liaison=dict()
    liaison2=dict()
    dic=dict()
    index_var_sup=dict()
    for i in range(len(b)):
        list_var.append("X{}".format(len(list_var)))
        var_artificielle.append("X{}".format(len(list_var)-1))
        index_var_sup[var_artificielle[-1]]=i
        if b[i]<0:
            b[i]=-b[i]
            index_const_neg.append(i)
            list_var.append("X{}".format(len(list_var)))
            var_exces.append("X{}".format(len(list_var)-1))
            cpt+=1
            dic[var_exces[-1]]=var_artificielle[-1]
            index_var_sup[var_exces[-1]]=i
        else:
            var_exces.append(None)
        liaison[var_artificielle[-1]]=var_exces[-1]
        liaison2[var_exces[-1]]=var_artificielle[-1]
    var_base=[]
    for i in range(len(var_exces)):
        if var_exces[i]!=None:
            var_base.append(var_exces[i])
    i=0
    while len(var_base)!=len(b):
        var_base.append(var_artificielle[i])
        i+=1
    c=np.copy(A)
    mat=[]
    btemp=[]
    for j in range(len(b)):
        btemp.append(b[index_var_sup[var_base[j]]])
        tmp=[]
        for i in range(len(f),len(list_var)):
            if list_var[i]==var_base[j]:
                tmp.append(1)
            elif var_base[j] in dic.keys():
                if liaison2[var_base[j]]==list_var[i]:
                    tmp.append(-1)
                else:
                    tmp.append(0)
            
            else:
                tmp.append(0)
        tmp=np.array(tmp)
        indice=index_var_sup[var_base[j]]
        if indice in index_const_neg:
            c[indice]=-c[indice]
        mat.append(np.concatenate((c[indice],tmp)))
    b=np.array(btemp)
    c=np.array(mat)
    return b,c,var_base,list_var,var_exces,var_artificielle

def init_Z(list_var,var_exces,f,b,var_base,c):
    Z=[]
    objec=0
    for i in range(len(list_var)):
        if len(var_exces)==0:
            Z=np.concatenate((f,np.zeros(len(list_var)-len(f))))
        else:
            somme=0
            for j in range(len(b)):
                if var_base[j] in var_exces:
                    if i ==0:
                        objec+=-b[j]
                    if i<len(f):
                        somme-=c[j][i]
                        
                    else:
                        if c[j][i]<0:
                            somme-=c[j][i]
            Z.append(somme)
    Z=np.array(Z)
    return Z,objec

def simplex_p1(itera,Z,c,b,var_base,list_var,VE,objec,verbose):
    """phase 1"""
    while True:
        index_entrant=find_max(Z)
        if index_entrant==None:
            break
        index_sortant=find_var_sortante(c,b,index_entrant)
        if index_sortant==None:
            break
        var_base[index_sortant]=list_var[index_entrant]
        
        c,Z,b,objec=pivot(list_var[index_entrant],var_base[index_sortant],c,var_base,list_var,VE,Z,b,objec)
        itera+=1
        if verbose:
            print("iteration={}".format(itera))
    return c,Z,b,objec

def transi_p1_p2(f,c,Z,list_var,var_exces,VE,objec,var_base,b):
    """transition between phase 1 and 2"""
    c=c.T
    tmp=[]
    cpt=0
    p=len(c)
    for i in range(p):
        if round(Z[i])!=1 and list_var[i] not in var_exces:
            cpt+=1
            tmp.append(c[i])
        else:
            VE.pop(VE.index(list_var[i]))
            list_var[i]=None
    c=np.array(tmp).T
    valeur_variable=dict()
    i=0
    list_var2=[]
    for a in list_var:
        if a!=None:
            list_var2.append(a)
            if i<len(f):
                valeur_variable[a]=f[i]
            else:
                valeur_variable[a]=0
            i+=1
    Z=np.zeros(cpt)
    i=0
    for a in list_var:
        if a!=None:
            somme=-valeur_variable[a]
            for j in range(len(var_base)):
                somme+=valeur_variable[var_base[j]]*c[j][list_var2.index(a)]
            Z[i]=somme
            i+=1
    for j in range(len(b)):
        objec+=valeur_variable[var_base[j]]*b[j]

    return var_base,list_var2,c,Z,objec

def simplex_p2(itera,Z,c,b,var_base,list_var2,VE,objec,tableau_objectif,verbose):
    """phase 2"""
    itera2=1
    while True:
        index_entrant=find_max(Z)
        if index_entrant==None:
            break
        index_sortant=find_var_sortante(c,b,index_entrant)
        if index_sortant==None:
            break        
        var_base[index_sortant]=list_var2[index_entrant]
        c,Z,b,objec=pivot(list_var2[index_entrant],var_base[index_sortant],c,var_base,list_var2,VE,Z,b,objec)
        tableau_objectif.append(-objec)
        itera+=1
        itera2+=1
        if verbose:
            print("iteration={}".format(itera))
    return Z,c,b,objec,tableau_objectif,itera2

def simplex(f,A,b, verbose=True):
    
    """minimise fx s.c Ax<=b
    Ce programme ne vérifie pas que le probleme a une solution.
    param:
    :verbose: Affiche l'évolution de la fonction objective au cours des itérations et le numéro de l'itération en cours 
    """
    f=-f

    b,c,var_base,list_var,var_exces,var_artificielle=init_simplex(f,A,b)
    Z,objec=init_Z(list_var,var_exces,f,b,var_base,c)
    tableau_objectif=[]
    VE=var_artificielle+var_exces 
    itera=0
    c,Z,b,objec=simplex_p1(itera,Z,c,b,var_base,list_var,VE,objec,verbose)


    if abs(objec) >10e-3:
        print("pas de solution")
        return#no solution

    var_base,list_var2,c,Z,objec=transi_p1_p2(f,c,Z,list_var,var_exces,VE,objec,var_base,b)
    tableau_objectif.append(-objec)
    itera+=1
    Z,c,b,objec,tableau_objectif,itera2=simplex_p2(itera,Z,c,b,var_base,list_var2,VE,objec,tableau_objectif,verbose)



    x=[]
    for i in range(len(f)):
        if list_var2[i] in var_base:
            x.append(b[var_base.index(list_var2[i])])
        else:
            x.append(0)
    if verbose:
        plt.plot(range(itera2),tableau_objectif)
        plt.xlabel("Iteration")
        plt.ylabel("Valeur de la fonction objectif")
        plt.title("Evolution de la valeur de la fonction objectif en fonction des itérations")
        #plt.legend()
        plt.show()
    
    return x
    
def main():
    #2 petits tests
    print("test 1")
    b = np.array([5,9]).T
    C = np.array([[-2,-3],[-4,-12],[-8,-6]])
    d = np.array([-90,-240,-240]).T
    res=simplex(b,C,d)#sol = x1=30,x2=10
    
    print("res")
    print(res)
    print("test 2")
    b = np.array([5,9]).T
    C = np.array([[-2,-3],[-4,-12],[-8,-6]])
    d = np.array([90,-240,-240]).T
    res=simplex(b,C,d)#sol = x1=20,x2=40/3
    print("res")
    print(res)
    
if __name__=="__main__":    
    main()
