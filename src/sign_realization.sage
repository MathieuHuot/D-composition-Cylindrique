#(***************************************************************************************)
#(*                                                                                     *)
#(*                                                                                     *)
#(*                              GARNIER Remy                                           *)
#(*                              HUOT Mathieu                                           *)
#(*                    Licence 3 : stage de Mathématiques                               *)
#(*                           Version Q[X1,...Xn]                                       *)
#(*                            Sign Realization                                         *)
#(*                                                                                     *)
#(*                                                                                     *)
#(***************************************************************************************)

import numpy #for tensorial product

#1°)Index conversion
#~~~~~~~~~~~~~~~~~~~

#INPUT : elem  'a
#        liste 'a list
#OUTPUT: i     integer : position of elem in liste and 0 if not found
def ConvertUplet(elem,liste):
    l=len(liste)
    for i in range(0,l):
        if liste[i] == elem :  
            return i
    return 0

#2°)Support et extraction matricielle
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#INPUT : M 'a matrix
#OUTPUT: E 'a list list : passe d'une matrice colonne a la liste de ses coefficients
def MatriceEnListe(M):
    l=M.nrows()
    E=[]
    for i in range(0,l):
        E=E+[M[i]]
    return E

#INPUT : M    'a matrix
#        ind  'b list : liste d'uplets
#OUTPUT: List 'b list : support de M matrice colonne en liste d'uplets
def Supp (M,ind):
    l=M.nrows()
    List=[]
    for i in range(0,l):
        if M[i]!=0:
            List=List+[ind[i]]
    return List

#INPUT : M        'a matrix
#        A        'b list : ensemble d'uplets
#        B        'b list : ensemble d'uplets
#        Lignes   'b list : indice les lignes de M par uplet
#        Colonnes 'b list : indice les colonnes de M par uplet
#OUTPUT: cMe      'a matrix : matrice extraite de M selon A et B
def Extract (M,A,B,Lignes,Colonnes):
    m=len(A)
    n=len(B)
    matrixe=MatrixSpace(QQ,m,n)
    cMe=[0 for i in range(0,m*n)]
    ic=0
    for UpletLigne in A:
        jc=0
        for UpletColonne in B:
            i=ConvertUplet(UpletLigne,Lignes)
            j=ConvertUplet(UpletColonne,Colonnes)
            cMe[ic*n+jc]=M[i][j]
            jc=jc+1
        ic=ic+1
    return matrixe(cMe)

#Extrait une suite L d'indices de M dont les lignes sont n premieres lignes 
#indépendantes de M
#INPUT : M        'a matrix
#        A        'b list : ensemble d'uplets
#        Sigma    'b list : ensemble d'uplets
#        n        integer
#        Lignes   'b list        : indice les lignes de M par uplet
#        Colonnes 'b list        : indice les colonnes de M par uplet
#OUTPUT: L        ('a list) list : n premières lignes de M indépendantes
def Extractlibre(M,A,Sigma,n,Lignes,Colonnes):
    L=[]
    i=0
    for ligne in range(0,len(A)):
        Ltest=L+[A[ligne]]
        itest=i+1
        Mtest=Extract(M,Ltest,Sigma,Lignes,Colonnes)
        if Mtest.rank()==itest:
            L=Ltest
            i=itest
        if i==n:
            return L
    return L 

#3°)Produit tensoriel:
#~~~~~~~~~~~~~~~~~~~~~

#INPUT : M1  'a matrix de taille m1*n1
#        m1  integer
#        n1  integer
#        M1  'a matrix de taille m2*n2
#        m2  integer
#        n2  integer
#OUTPUT: res 'a matrix de taille m1m2*n1n2 : produit tensoriel de M1 par M2
def Tproduit(M1,m1,n1,M2,m2,n2):
    return (MatrixSpace(QQ,m1*m2,n1*n2))((numpy.kron(M1,M2)).tolist())

#4°)Produits cartésiens 
#~~~~~~~~~~~~~~~~~~~~~~~

#INPUT : A 'a list 
#        B 'b list
#OUTPUT: E ('a * 'b) list : produit cartésien de A par B
def ProduitCart(A,B):
    E=[]
    i=0
    for ela in A:
        for elb in B:
            E=E+[ela+elb]    
    return E

#5°) Algorithme proprement dit:
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#INPUT : l     integer 
#        T     (int * Q[X1,...,Xl] * int) list : système triangulaire
#        P     Q[X1,...,Xl-1][Xl]
#        p     integer : degré de P
#        Q     Q[X_1,...,X_l] list
#        deg   int list
#OUTPUT: Sigma (int list) list 
#        nb    (int list) list
def SignRealization(l,T,P,p,Q,deg):
#5.1: Calculs préliminaires
#On définit d'abord les espaces de matrices initiaux, ainsi que les matrices de départ 
    MS=MatrixSpace(QQ,3,3)
    Mcol=MatrixSpace(QQ,3,1)
    M1=MS([1,1,1,-1,0,1,1,0,1])
    invM1=M1^(-1)
    M=M1
    #Initialisation des différents TaQ (sous forme de liste et de matrice (TaQM))    
    TaQM=Mcol([0,0,0]) 
    TaQ=Zero(3)
    m=len(Q)#nombre de polynomes Q
    dP=TdA[l](diff(TdA[l](P),TdV[l-1]))#Dérivée de P
    #dP=Primitif(l,dP)#On le rend primitif
#5.2: Cas de Base   
    for e in [0,1,2]:
        Ro=dP*(Q[m-1])^e
        r=p-1+deg[m-1]*e
        R,r=IntRem(l,T,Ro,r,P,p)
        #R=Primitif(l,R) #On le rend primitif
        TaQ[e]=PmVPol(l,T,P,p,R,r) #Calcul des premiers TaQ
    TaQM=Mcol(TaQ)
    nbM=invM1*TaQM
    if nbM==0:
        return [] #P est de signe constant
    Sigma=Supp(nbM,[[0],[1],[2]])
    lep=len(Sigma)
    if lep==1:
        A=[[0]]
    elif lep==2:
        A=[[0],[1]]
    else:
        A=[[0],[1],[2]]
    M=Extract(M1,A,Sigma,[[0],[1],[2]],[[0],[1],[2]])
    nb=Extract(nbM,Sigma,[[0]],[[0],[1],[2]],[[0]])
    AncienSigma=Sigma #Servent à indexer lors de l'extraction de lignes libres de A
    AncienA=A
#5.3: Boucle principale
    for i in range(m-1,0,-1) :
        extSigma=ProduitCart([[0],[1],[2]],Sigma) 
#A priori le fait de numéroter par {0,1,2}, plutot que {-1,0,1} ne change rien 
#et permet d'utiliser les memes fonctions que pour A .
        extA=ProduitCart([[0],[1],[2]],A)
        extM=Tproduit(M1,3,3,M,M.nrows(),M.ncols())
        Mcoli=MatrixSpace(QQ,len(extA),1)
        TaQi=Zero(len(extA))
        for Uplet in extA :
            R=dP
            r=p-1
            for j in range(i-1,m): 
                R=R*(Q[j])**(Uplet[j-i+1])
                r=r+deg[j]*(Uplet[j-i+1])
            R,r=IntRem(l,T,R,r,P,p)
            #R=Primitif(l,R) #On le rend primitif
            TaQi[ConvertUplet(Uplet,extA)]=PmVPol(l,T,P,p,R,r)
        TaQM=Mcoli(TaQi)
        nbM=((extM)^(-1))*TaQM
        Sigma1=Sigma
        Sigma=Supp(nbM,extSigma)
        Sigma2=[]
        Sigma3=[]
        for sigma in Sigma1:
            compteur=0
            for i in range(0,3):
                if [i]+sigma in Sigma:
                    compteur += 1
            if compteur>2:
                Sigma3=[sigma]+Sigma3
            if compteur>1:
                Sigma2=[sigma]+Sigma2
        s2=len(Sigma2)
        s3=len(Sigma3)
        A1=A
        A2=Extractlibre(M,A,Sigma2,s2,AncienA,AncienSigma)
        A3=Extractlibre(M,A2,Sigma3,s3,AncienA,AncienSigma)
        A=(ProduitCart([[0]],A1))+(ProduitCart([[1]],A2))+(ProduitCart([[2]],A3))
        nb=Extract(nbM,Sigma,[Zero(m+1-i)],extSigma,[Zero(m+1-i)])        
        M=Extract(extM,A,Sigma,extA,extSigma)
        AncienSigma=Sigma
        AncienA=A
    nb=MatriceEnListe(nb)#Eh oui, nb est une matrice, et non une liste
    Sigma=map(DecrL,Sigma)#On retourne à un codage dans {-1,0,1}
    return Sigma,nb 
