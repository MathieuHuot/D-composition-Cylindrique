#(***************************************************************************************)
#(*                                                                                     *)
#(*                                                                                     *)
#(*                              HUOT Mathieu                                           *)
#(*                              GARNIER Remy                                           *)
#(*                    Licence 3 : stage de Mathématiques                               *)
#(*                           Version Q[X1,...Xn]                                       *)
#(*                            Sign Realization                                         *)
#(*                                                                                     *)
#(*                                                                                     *)
#(***************************************************************************************)

#1°)Conversion des indices
#~~~~~~~~~~~~~~~~~~~~~~~~~

#Renvoie la position d'un élement dans une liste et 0 si non trouvé
def ConvertUplet(elem,liste):
    l=len(liste)
    for i in range(0,l):
        if liste[i] == elem :  
            return i
    return 0
#'a * 'a list -> int

#2°)Support et extraction matricielle
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#Passer d'une matrice colonne à la liste de ses coefficients:
def MatriceEnListe(M):
    l=M.nrows()
    E=[]
    for i in range(0,l):
        E=E+[M[i]]
    return E
#'a Matrix -> ('a list) list

#Renvoie le support d'une matrice colonne sous forme d'une liste d'uplets. 
def Supp (M,indices):
    l=M.nrows()
    List=[]
    for i in range(0,l):
        if M[i]!=0:
            List=List+[indices[i]]
    return List
#'a Matrix * 'b list -> 'b list    

#Construit une sous matrice de M à partir d'un ensemble de lignes A et d'un ensemble 
#de colonne B. La fonction prend également en paramètre deux ensembles d'indices 
#Lignes et Colonnes contenant respectivement A et B et indiciant M.   
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
#'a Matrix * 'b list * 'b list * 'b list * 'b list -> 'a Matrix

#Extrait une suite L d'indices de M dont les lignes sont n premieres lignes 
#indépendantes de M
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
# 'a Matrix *   -> ('a list) list

#3°)Produit tensoriel:
#~~~~~~~~~~~~~~~~~~~~~

#Cet fonction réalise le produit tensoriel de 2 matrices M1 et M2, de taille respectives
# (m1,n1) et (m2,n2)
def Tproduit(M1,m1,n1,M2,m2,n2):
    gromatrix=MatrixSpace(QQ,m1*m2,n1*n2)
    nMe=[0 for i in range(0,m1*m2*n1*n2)]
    for i1 in range(0,m1):
        for i2 in range(0,m2):
            for j1 in range(0,n1):
                for j2 in range(0,n2):
                    nMe[(i1*m2+i2)*n1*n2+j1*n2+j2]=(M1[i1][j1])*(M2[i2][j2])
    return gromatrix(nMe)
#'a Matrix * int * int * 'a Matrix * int * int -> 'a Matrix

#4°)Produits cartésiens 
#~~~~~~~~~~~~~~~~~~~~~~~

#Effectue en terme de liste le produit cartésien de 2 ensembles A et B
def ProduitCart(A,B):
    E=[]
    i=0
    for ela in A:
        for elb in B:
            E=E+[ela+elb]    
    return E
#'a list * 'b list  -> ('a * 'b) list

#5°) Algorithme proprement dit:
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#P est dans D[X], D=Q[X_1,...X_l-1], X=X_l
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
    dP=diff(P+0*TdV[l-1],TdV[l-1])#Dérivée de P
    dP=Primitif(l,dP)#On le rend primitif
#5.2: Cas de Base   
    for e in [0,1,2]:
        Ro=dP*(Q[m-1])^e
        r=p-1+deg[m-1]*e
        R,r=IntRem(l,T,Ro,r,P,p)
        R=Primitif(l,R) #On le rend primitif
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
            R=Primitif(l,R) #On le rend primitif
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
#int * (int * Q[X_1,...,X_l] * int) list *  Q[X_1,...,X_l] * int *  Q[X_1,...,X_l] list *
#int list -> (int list) list * (int list) list
