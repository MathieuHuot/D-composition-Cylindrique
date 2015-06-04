#(***************************************************************************************)
#(*                                                                                     *)
#(*                                                                                     *)
#(*                              HUOT Mathieu                                           *)
#(*                              GARNIER Remy                                           *)
#(*                    Licence 3 : stage de Mathématiques                               *)
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
    if a[k+1]==1 and a[k]<b[k]:
        return -1
    if a[k+1]==-1 and a[k]>b[k]:
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
        res= res + [[[i+1,L[i],P]]] #on code un triplet par une liste
    return res
#COMPLEXITY : O(len(L))

#Agrandit chaque élément Sli[j] de Sli par (r,SLL[j],Q) où r désigne le numéro de 
#la racine si SLL[j] code une racine de Q et -1 sinon
#Version où on concatène à gauche
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
#Version où on concatène à droite
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
        res=res + [SLj[k]+[[r,SLL[k],Q]]]
    return res
#((int * (int list) * Q[X1,...Xl]) list) list * (int list) list * Q[X1,...Xl] *
#((int * (int list) * Q[X1,...Xl]) list) list * int 
# -> ((int * (int list) * Q[X1,...Xl]) list) list

#La comparaison utilisée pour le tri nécessaire dans OrderedMerge . On fait appel à 
#PetitThom2 sur deux P-encoding dont l'un au moins est sur une racine de P, ce qui 
#garantit la validité de la comparaison totale de la liste en ne comparant que deux 
#éléments, et on trouve un tel P car on a sauvegardé son indice en début de liste
def CompP(L1,L2):
    i=L1[0]
    return PetitThom2(L1[i][1],L2[i][1])
#(int * int list * Q[X1,...,Xl]) list * (int * int list * Q[X1,...,Xl]) list -> {-1,0,1}

#Effectue le tri optimisé  de deux listes triées pour la fonction OrderedMerge
def TriP(L1,L2):
    n=len(L1)
    m=len(L2)
    res=[]
    i=0
    j=0
    while(i<n and j<m):
        k=CompP(L1[i],L2[j])
        if k==-1:
            res=res+[L1[i]]
            i+=1
        elif k==1:
            res=res+[L2[j]]
            j+=1
        else:
            res=res+[L1[i]]
            i+=1
            j+=1
    for k in range(i,n):
        res=res+[L1[k]]
    for k in range(j,m):
        res=res+[L2[k]]
    return res        
#(int * int list * Q[X1,...,Xl]) list list * (int * int list * Q[X1,...,Xl]) list list ->
#(int * int list * Q[X1,...,Xl]) list list

#Effectue la fusion d'une liste de listes en une liste, qui est triée et dont
#on a enlevé les doublons : c'est à dire trie les racines des différents polynomes à
#l'aide de leur Thom-Encoding
def OrderedMerge(SL): 
    n=len(SL)
    mid=n//2
    if n==0:
        return []
    elif n==1: #############rajouter le cas if n==0 renvoyer [] fait aussi planter
        return SL[0]
    elif n==2:
        return TriP(SL[0],SL[1])
  
    else:
        L1=[SL[i] for i in range(mid)]
        L1=OrderedMerge(L1)
        L2=[SL[i] for i in range(mid,n)]
        L2=OrderedMerge(L2)
        return TriP(L1,L2)
#(int * int list * Q[X1,...,Xl]) list list list -> 
#(int * int list * Q[X1,...,Xl]) list list

#Renvoie le polynome sans ses coefficients nuls lorsque précisé en alpha_1,...,alpha_l
def Normalize(l,T,P):
    p = Degree(l,T,P)
    Res = 0
    for j in range(0,p+1):
        Res = Res + P[j]*TdV[l-1]**j #La variable principale de P est X_l
    Res = TdA[l](Res*(1/A(Res).content()))    #On rend le polynome primitif
    return Res,p
#int * (int * Q[X1,...Xl] * int) list * Q[X1,...Xl] -> Q[X1,...Xl]

#Effectue la partition de la ligne réelle selon les polynomes de PP à l variables
#dans le système triangulaire T de niveau l-1
def LinePartition(PP2,l,T):
    lon2=len(PP2) #PP2 a des polynomes à l variables
    Normed2=[Normalize(l,T,PP2[i]) for i in range(lon2)] #Préparcours pour calculer
    Normed=[] #les polynomes normalisés avec racines qui vont aller dans Normed
    RootCodePi=[] #et leur Rootcoding dans RootCodePi
    List=[]
    for i in range(long2):
        if pi>0:
            Pi,pi=Normed2[i]
            List=List+ [(l,T,)]
    
    Output=RootPar(List) # Parallelisation des rootcodings
    
    for i in range(len(List)):
        Root=Output[i][2]
        if Root!=[]: #On ne garde que les polynomes qui ont des racines
            Normed=Normed+[Normed2[i]]
            RootCodePi=RootCodePi+[Root]



    lon=len(Normed)#On place tous les RootCoding des polynomes normalisés dans une matrice
    SLL=[[0 for j in range(lon)] for i in range(lon)]
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
#Q[X1,...,Xl] list * int * (int * Q[X1,...,Xl] * int)  list -> 
#(int * int list * Q[X1,...,Xl]) list list * (Q[X1,...,Xl] * int) list
