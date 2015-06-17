#(***************************************************************************************)
#(*                                                                                     *)
#(*                                                                                     *)
#(*                              HUOT Mathieu                                           *)
#(*                              GARNIER Remy                                           *)
#(*                    Licence 3 : stage de Mathématiques                               *)
#(*                           Version Q[X1,...Xn]                                       *)
#(*                    Completing the real line partition                               *)
#(*                                                                                     *)
#(*                                                                                     *)
#(***************************************************************************************)

#INPUT : lis  ((int * (int list) * Q[X_1,...,X_l]) list) list 
#        code (int list) list 
#        Q    Q[X_1,...,X_l] 
#OUTPUT: res  ((int * (int list) * Q[X_1,...,X_l]) list) list
#NOTE : Version plus basique qui est un EnlargeWithGauche sans recherche de racine
def EnlargeWithCompleting(lis,code,Q):
    n=len(lis)
    res=[]
    for i in range(n):
        res=res+[[[-1,code[i],Q]]+lis[i]]            
    return res
#COMPLEXITY : O(len(lis))

#INPUT : shortL (int * int list * Q[X_1,...,X_l]) list : liste de codage de racine
#        v      int list : Thom Encoding de P pour une racine \alpha_i
#        oldv   int list : Thom Encoding de oldP pour une racine \alpha_j < \alpha_i
#        P      Q[X_1,...,X_l]
#        oldP   Q[X_1,...,X_l]
#OUTPUT: F      int * (int * int list * Q[X_1,...,X_l]) list 
#               : racine dans shortL (vp,P) tel que vp<v voldP>oldv pour le Thom Encoding
def Find(shortL,v,oldv,P,oldP):
    n=len(shortL)
    r1=shortL[0]
    m=len(r1)
    indiceP=-1
    indiceOP=-1
    for i in range(m): #On trouve les indices i et j tels que P=Pi, oldP=Pj
        if P==r1[i][2]:
            indiceP=i
            break
    for i in range(m):
        if oldP==r1[i][2]:
            indiceOP=i
            break
    for i in range(n): #On cherche la racine voulue F parmi les racines L[i] dans L
        F=shortL[i]
        vP=F[indiceP][1]
        voldP=F[indiceOP][1]
        if PetitThom(vP,v) and PetitThom(oldv,voldP):
            Ind=len(F) #On utilise la meme technique qu'avant : le 1er élément est 
            return [Ind]+F #l'indice de F où "r est défini"
#COMPLEXITY : O(len(shortL) + len(shortL[0])) 

#INPUT : l  integer
#        T  (int * Q[X1,...,Xl] * int) list 
#        PP (Q[X1,...,Xl] * int) list 
#        i  integer
#        j  integer
#OUTPUT: shortL (int * int list * Q[X1,...,Xl]) list list : R-rootcoding des racines de (PiPj)'
def PreCalculCompleting(l,T,PP,i,j):
    m=len(PP)
    P,p=PP[i]
    Q,q=PP[j]
    pol=diff(P*Q,TdV[l-1])
    pol=Primitif(l,pol) #On le rend primitif
    pol,p=Normalize(l,T,pol)
    SLL=RootCoding(l,T,pol,p,pol,p)
    shortL=Singleton(SLL,pol) 
    for k in range(m-1,-1,-1):
        R,r=PP[k]
        SLL=RootCoding(l,T,pol,p,R,r)
        shortL=EnlargeWithCompleting(shortL,SLL,R)
    return shortL
#COMPLEXITY : O(2EXP)

#INPUT : l  integer
#        T  (int * Q[X1,...,Xl] * int) list 
#        PP (Q[X1,...,Xl] * int) list 
#        i  integer
#OUTPUT: shortL (int * int list * Q[X1,...,Xl]) list list : R-rootcoding des racines de Pi'
def PreCalculCompleting2(l,T,PP,i):
    m=len(PP)
    P,p=PP[i] #Le polynome Pi et son degré
    P=diff(P,TdV[l-1]) #Pi'
    p-=1 #Et son degré qui vaut le degré de Pi -1
    P=Primitif(l,P) #On le rend primitif
    SLL=RootCoding(l,T,P,p,P,p)
    shortL=Singleton(SLL,P) 
    for j in range(m-1,-1,-1):
        Q,q=PP[j]
        SLL=RootCoding(l,T,P,p,Q,q)
        shortL=EnlargeWithCompleting(shortL,SLL,Q)
    return shortL
#COMPLEXITY : O(2EXP)

#INPUT : l    integer
#        T    (int * Q[X1,...,Xl] * int) list 
#        PP   (Q[X1,...,Xl] * int) list 
#        oldP Q[X1,...,Xl] 
#        P    Q[X1,...,Xl] 
#OUTPUT: shortL (int * int list * Q[X1,...,Xl]) list list : R-rootcoding des racines de (PoldP)'
def CalculCompleting(l,T,PP,oldP,P):
    m=len(PP)
    pol=diff(P*oldP,TdV[l-1])
    pol=Primitif(l,pol) #On le rend primitif
    pol,p=Normalize(l,T,pol)
    SLL=RootCoding(l,T,pol,p,pol,p)
    shortL=Singleton(SLL,pol) 
    for i in range(m-1,-1,-1):
        Q,q=PP[i]
        SLL=RootCoding(l,T,pol,p,Q,q)
        shortL=EnlargeWithCompleting(shortL,SLL,Q)
    return shortL
#COMPLEXITY : O(2EXP)

#INPUT : l    integer
#        T    (int * Q[X1,...Xl] * int) list 
#        L    (int * int list * Q[X1,...,Xl]) list list 
#        PP   (Q[X1,...,Xl] * int) list 
#OUTPUT: newL (int * int list * Q[X1,...,Xl]) list list
#NOTE : Complète la ligne réelle de niveau l partitionnée par line partition en ajoutant 
#une racine pour représenter chaque cellule qui n'est pas un singleton (les singletons 
#étant les racines des polynomes obtenus durant la phase d'élimination)
def Completing(l,T,L,PP):
    n=len(L) 
    m=len(PP)
#Il y a n racines donc 2n+1 représentants des cellules qu'on stockera dans newL
    newL=Zero(2*n+1) 
#On va précalculer tous les rootcoding des (PiPj)' où l'on a : L[i]=racine de (Pi ou Pj)
#et L[i+1]=racine de (Pi ou Pj) et cela pour au moins deux "i" distincts : cela évite
#par exemple le calcul de (PiPj)' puis (PjPi)' ou encore de (PiPi)'.
#On stocke dans Mtemp les indices i,j des PiPj à précalculer et dans M les résultats des
#calculs.
    M=[["Vide" for i in range(m)] for j in range(m)] #pour stocker les Root-Coding
    Mtemp=[[0 for i in range(m)] for j in range(m)] #pour compter si un précalcul est rentable
    
    for i in range(n-1):
        e=L[i]
        f=L[i+1]
        Ind1=e[0]-1 #P_Ind1(e)=0
        Ind2=f[0]-1 #P_Ind2(f)=0
        Mtemp[Ind1][Ind2]+=1
    if PAR_PCC:
        MPar=[] #pour pouvoir paralléliser les précalculs
        for i in range(m): 
            for j in range(i):
                Mtemp[i][j]+=Mtemp[j][i]
                if Mtemp[i][j]>1: #On a alors intéret à précalculer car le cas ci-dessus 
                    MPar+=[(l,T,PP,i,j)] #arrive au moins deux fois
        ResPar=list(ParPreCalculCompleting(MPar))
        for k in range(len(ResPar)):
            i,j=ResPar[k][0][0][3],ResPar[k][0][0][4]
            M[i][j]=ResPar[k][1]
    else:
        for i in range(m): 
            for j in range(i):
                Mtemp[i][j]+=Mtemp[j][i]
                if Mtemp[i][j]>1:  
                    M[i][j]=PreCalculCompleting(l,T,PP,i,j)
    
    if PAR_PCC2:    
        MPar=[]
        for i in range(m):
            if Mtemp[i][i]>0: #ie on veut éviter de calculer (PiPi)'
                MPar+=[(l,T,PP,i)]
        ResPar=list(ParPreCalculCompleting2(MPar))
        for k in range(len(ResPar)):
            i=ResPar[k][0][0][3]
            M[i][i]=ResPar[k][1]
    else:
        for i in range(m):
            if Mtemp[i][i]>0:
                M[i][i]=PreCalculCompleting2(l,T,PP,i)
    
    for i in range(n): #On recopie ce qui ne change pas : les cellules qui sont des
        newL[2*i+1]=L[i] #singletons contenant les racines des polynomes de PP

    E=L[0] #Cas de la cellule où il y a -l'infini
    r,v,P=E[E[0]] #E[0] donne l'indice du polynome dont E code une racine
    R,r=Normalize(l,T,Decale(l,P,1))
    SLL=RootCoding(l,T,R,r,R,r)
    newL[0]=Singleton(SLL,R)
    
    if ROOT_PAR:
        MPar=[]
        for j in range(m-1,-1,-1):
            Q,q=PP[j]
            MPar+=[(l,T,R,r,Q,q,j)]
        ResPar=list(RootPar(MPar))
        for j in range(m-1,-1,-1):
            for k in range(len(ResPar)):
                if ResPar[k][1][1]==j:
                    Q,q=PP[j]
                    newL[0]=EnlargeWithCompleting(newL[0],ResPar[k][1][0],Q)
    else:
        for j in range(m-1,-1,-1):
        Q,q=PP[j]
        SLL=RootCoding(l,T,R,r,Q,q)
        newL[0]=EnlargeWithCompleting(newL[0],SLL,Q)
                    
    newL[0]=[len(newL[0][0])]+newL[0][0] #C'est la 1ère racine qui nous intéresse
    oldv=v
    oldP=P

    for i in range(1,n): #Cas général de la cellule entre deux racines
        E=L[i]
        r,v,P=E[E[0]]
        Ind1=max(L[i-1][0],L[i][0])-1
        Ind2=min(L[i-1][0],L[i][0])-1
        if Ind1==Ind2: #On regarde si on a déjà précalculé le résultat
            newL[2*i]=M[Ind1][Ind2]
        elif Mtemp[Ind1][Ind2]>1: 
            newL[2*i]=M[Ind1][Ind2]
        else: #Si non, on le calcule maintenant
            newL[2*i]=CalculCompleting(l,T,PP,oldP,P)
        #On stocke la bonne racine ie celle entre L[i] et L[i+1]
        newL[2*i]=Find(newL[2*i],v,oldv,P,oldP)
        oldv=v
        oldP=P

    E=L[n-1] #Cas de la cellule où il y a +l'infini
    r,v,P= E[E[0]]
    R,r=Normalize(l,T,Decale(l,P,-1))
    SLL=RootCoding(l,T,R,r,R,r)
    newL[2*n]=Singleton(SLL,R)
    
    if ROOT_PAR:
        MPar=[]
        for j in range(m-1,-1,-1):
            Q,q=PP[j]
            MPar+=[(l,T,R,r,Q,q,j)]
        ResPar=list(RootPar(MPar))
        for j in range(m-1,-1,-1):
            for k in range(len(ResPar)):
                if ResPar[k][1][1]==j:
                    Q,q=PP[j]
                    newL[2*n]=EnlargeWithCompleting(newL[2*n],ResPar[k][1][0],Q)
    else:
        for j in range(m-1,-1,-1):
            Q,q=PP[j]
            SLL=RootCoding(l,T,R,r,Q,q)
            newL[2*n]=EnlargeWithCompleting(newL[2*n],SLL,Q)  
            
    #C'est la plus grande racine qui nous intéresse :
    newL[2*n]=newL[2*n][len(newL[2*n])-1]
    newL[2*n]=[len(newL[2*n])]+newL[2*n] #On rajoute encore le "i" initial indiquant où chercher le
    return newL                          # "r" défini
#COMPLEXITY : O(2EXP)
