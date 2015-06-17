#(***************************************************************************************)
#(*                                                                                     *)
#(*                                                                                     *)
#(*                              HUOT Mathieu                                           *)
#(*                              GARNIER Remy                                           *)
#(*                    Licence 3 : stage de MathÃ©matiques                               *)
#(*                           Version Q[X1,...Xn]                                       *)
#(*                              Line Partition                                         *)
#(*                                                                                     *)
#(*                                                                                     *)
#(***************************************************************************************)

#INPUT : a   int list
#        b   int list
#OUTPUT: res {-1,0,1} : 1 si a<b, 0 si a=b, -1 si b>a pour l'ordre du Thom Encoding
def PetitThom2(a,b):
    if a==b:
        return 0
    k=len(a)-1
    while a[k]==b[k]:
        k-=1
    if (a[k+1]==1 and a[k]<b[k]) or (a[k+1]==-1 and a[k]>b[k]):
        return -1
    return 1
#COMPLEXITY : WORSE : O(len(a)) AVERAGE : O(1)

#INPUT : L   (int list) list 
#        P   Q[X1,...Xl] 
#OUTPUT: res ((int * (int list) * Q[X1,...Xl]) list) list
#NOTE  : Transforme une liste d'objets en une liste de singletons contenant 
#        ces objets et ajoute le numéro de la racine de P
def Singleton(L,P):
    res=[]
    l=len(L)
    for i in range(0,l):
        res+=[[[i+1,L[i],P]]] #on code un triplet par une liste
    return res
#COMPLEXITY : O(len(L))

#INPUT : SLj     ((int * (int list) * Q[X1,...Xl]) list) list
#        SLL     (int list) list
#        NormedQ Q[X1,...Xl]
#        SLi     ((int * (int list) * Q[X1,...Xl]) list) list
#        i       integer
#OUTPUT: res     ((int * (int list) * Q[X1,...Xl]) list) list
#NOTE  : agrandit chaque élément Sli[j] de SLi par (r,SLL[j],Q) où r désigne le numéro de la racine
#        si SLL[j] code une racine de Q et -1 sinon.
#        version où on concatène à gauche
def EnlargeWithGauche(SLi,SLL,NormedQ,SLj,J):
    res=[]
    lon=len(SLL)
    lon2=len(SLj)
    lon3=len(SLi)
    Q,q=NormedQ
    for j in range(lon):
        r=-1 #Par défaut ce n'est pas une racine de Q
        for k in range(lon2): #On cherche si un Q-encoding est identique
            if SLL[j]==SLj[k][J][1]: #Dans ce cas c'est une racine de Q
                r=k+1
                break
        res=res + [[[r,SLL[j],Q]]+SLi[j]]
    for j in range(lon3,lon):
        res=res+[SLi[j]]
    return res
#COMPLEXITY : O(lon * lon2)    
    
#INPUT : SLj     ((int * (int list) * Q[X1,...Xl]) list) list
#        SLL     (int list) list
#        NormedQ Q[X1,...Xl]
#        SLi     ((int * (int list) * Q[X1,...Xl]) list) list
#        i       integer
#OUTPUT: res     ((int * (int list) * Q[X1,...Xl]) list) list
#NOTE  :  version où on concatène à droite
def EnlargeWithDroite(SLj,SLL,NormedQ,SLi,i):
    res=[]
    lon=len(SLL)
    lon2=len(SLi)
    Q,q=NormedQ
    for k in range(lon):
        r=-1
        for g in range(lon2):
            if SLL[k]==SLi[g][i][1]:
                r=g+1
                break
        res+=[SLj[k]+[[r,SLL[k],Q]]]
    return res
#COMPLEXITY : O(lon * lon2)   

#INPUT : L1 (int * int list * Q[X1,...,Xl]) list 
#        L2 (int * int list * Q[X1,...,Xl]) list 
#OUTPUT: s  {-1,0,1} comparaison de L1 et L2 codant une racine chacun pour l'ordre du Thom Encoding
#NOTE  : on fait appel à PetitThom2 sur 2 P-encoding dont l'un au moins est sur une racine de P 
#        (l'indice au début est la pour ça), ce qui garantit la validité de la comparaison totale
#        de la liste en ne comparant que deux éléments
def CompP(L1,L2):
    i=L1[0]
    return PetitThom2(L1[i][1],L2[i][1])
#COMPLEXITY: AVERAGE O(1)

#INPUT : L1  (int * int list * Q[X1,...,Xl]) list list
#        L2  (int * int list * Q[X1,...,Xl]) list list
#OUTPUT: res (int * int list * Q[X1,...,Xl]) list list tri fusion optimisé de L1 et L2 
def TriP(L1,L2):
    n=len(L1)
    m=len(L2)
    res=[]
    i=0
    j=0
    while(i<n and j<m):
        k=CompP(L1[i],L2[j])
        if k==-1:
            res+=[L1[i]]
            i+=1
        elif k==1:
            res+=[L2[j]]
            j+=1
        else:
            res+=[L1[i]]
            i+=1
            j+=1
    for k in range(i,n):
        res+=[L1[k]]
    for k in range(j,m):
        res+=[L2[k]]
    return res        
#COMPLEXITY : O(n+m)

#INPUT : SL (int * int list * Q[X1,...,Xl]) list list list 
#OUTPUT: L  (int * int list * Q[X1,...,Xl]) list list : fusion des listes de SL avec TriP en enlevant
#           les doublons: i.e trie les racines des différents polynomes à l'aide de leur Thom-Encoding
def OrderedMerge(SL): 
    n=len(SL)
    mid=n//2
    if n==0:
        return []
    elif n==1: 
        return SL[0]
    elif n==2:
        return TriP(SL[0],SL[1])
    else:
        L1=OrderedMerge([SL[i] for i in range(mid)])
        L2=OrderedMerge([SL[i] for i in range(mid,n)])
        return TriP(L1,L2)
#COMPLEXITY : O( n * m * log(m) )

#INPUT : l   integer
#        T   (int * Q[X1,...,Xl] * int) list
#        P   Q[X1,...,Xl]
#OUTPUT: Res Q[X1,...,Xl-1][Xl] P sans ses coeffs nuls lorsque précisé en \alpha_1,...,\alpha_l-1
#        p   integer
def Normalize(l,T,P):
    p=Degree(l,T,P)
    Res=0
    for j in range(0,p+1):
        Res+=P[j]*TdV[l-1]**j #La variable principale de P est X_l
    Res=Primitif(l,Res)    #On rend le polynome primitif
    return Res,p
#COMPLEXITY: O(2EXP)

#INPUT : PP2    Q[X1,...,Xl] list 
#        l      integer 
#        T      (int * Q[X1,...,Xl] * int)  list : système triangulaire
#OUTPUT: SL     (int * int list * Q[X1,...,Xl]) list list 
#        Normed (Q[X1,...,Xl] * int) list : polynomes de PP2 normalisés
#NOTE  : SL représente la partition de la ligne réelle en étant une liste de codage de racines 
#        échantillon de chaque intervalle de la partition
def LinePartition(PP2,l,T):
    lon2=len(PP2) #PP2 a des polynomes à l variables
    
    Normed2=Zero(lon2)
    if NORM_PAR:
        Par=[] #liste à paralleliser
        for i in range(lon2):
            Par=Par+[(l,T,PP2[i],i)]
        Output=list(NormalizePar(Par))
        for i in range(lon2):
            j=Output[i][1][1]
            Normed2[j]=Output[i][1][0]
    else:
        Normed2=[Normalize(l,T,PP2[i]) for i in range(lon2)]
        
    Normed=[] #les polynomes normalisés avec racines qui vont aller dans Normed    
    if ROOT_PAR:
        RootCodePi=[] #et leur Rootcoding dans RootCodePi
        List=[]
        for i in range(lon2):
            Pi,pi=Normed2[i]
            if pi>0:
                List=List+ [(l,T,Pi,pi,Pi,pi,i)]
        lon=len(List)
        Output=list(RootPar(List)) # Parallelisation des rootcodings
        for i in range(lon):
            k=0
            for j in range(lon):
                if Output[j][1][1]==i:
                    k=j
                    break
            Root=Output[k][1][0]
            if Root!=[]: #On ne garde que les polynomes qui ont des racines
                Normed=Normed+[Normed2[i]]
                RootCodePi=RootCodePi+[Root]
    else:
        for i in range(lon2):
            Pi,pi=Normed2[i]
            if pi>0: #On ne veut pas des polynomes constants
                Root=RootCoding(l,T,Pi,pi,Pi,pi)           
                if Root!=[]: #On ne garde que les polynomes qui ont des racines
                    Normed=Normed+[Normed2[i]]
                    RootCodePi=RootCodePi+[Root]
            
    lon=len(Normed)#On place tous les RootCoding des polynomes normalisés dans une matrice
    SLL=[[0 for j in range(lon)] for i in range(lon)]
    
    if ROOT_PAR2:
        ListArg=[] #Liste des arguments pour la parallélisation
        for i in range(lon):
            for j in range(lon):
                if i==j: #On a déjà calculé et stocké le RootCoding de Pi sur ses racines
                    SLL[i][i]=RootCodePi[i]
                else:
                    Pi,pi=Normed[i]
                    Pj,pj=Normed[j]
                    ListArg=ListArg+[(l,T,Pi,pi,Pj,pj,i,j)]
        Output=list(RootPar2(ListArg))
        tal=len(Output)
        for i in range(lon):
            for j in range(lon):
                if i !=j:    
                    for k in range(tal):
                        if Output[k][1][1]==i and Output[k][1][2]==j:
                            SLL[i][j]=Output[k][1][0]
                            break
    else:
        for i in range(lon):
             for j in range(lon):
                if i==j: #On a déjà calculé et stocké le RootCoding de Pi sur ses racines
                     SLL[i][i]=RootCodePi[i]
                else:
                    Pi,pi=Normed[i]
                    Pj,pj=Normed[j]
                    SLL[i][j]=RootCoding(l,T,Pi,pi,Pj,pj)
    
    SL=[[] for i in range(lon)] 
    for i in range(lon):
        Pi,pi=Normed[i]
        SL[i]=Singleton(SLL[i][i],Pi) 
        for j in range(i-1,-1,-1):
            SL[i]=EnlargeWithGauche(SL[i],SLL[i][j],Normed[j],SL[j],j) 
            #On a besoin de SL[j] et j pour savoir si la racine est une racine de Q
        for j in range(i-1,-1,-1):
            SL[j]=EnlargeWithDroite(SL[j],SLL[j][i],Normed[i],SL[i],i)
    for i in range(lon):
        for j in range(len(SL[i])):
            SL[i][j]=[i+1]+SL[i][j] 
        #On précise devant chaque codage de racine le numéro du polynome qui l'engendre
    SL=OrderedMerge(SL)
    #for i in range(len(SL)): #On enlève l'information précédemment rajoutée
    #    SL[i]=[SL[i][j] for j in range(1,len(SL[i]))] 
    return SL,Normed
#COMPLEXITY : O(2EXP)
