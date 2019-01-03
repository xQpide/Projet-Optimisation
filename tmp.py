def init_simplex(f,A,b):
    tableau_objectif=[]
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
    return b,c,var_base,list_var,var_exces

def init_Z():
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
    return Z