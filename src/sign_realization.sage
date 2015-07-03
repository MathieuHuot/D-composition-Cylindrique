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

#2°)Support and matrix extraction
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#INPUT : M 'a matrix
#OUTPUT: E 'a list list : turns a column matrix into the list of its coefficients
def MatriceEnListe(M):
    l=M.nrows()
    E=[]
    for i in range(0,l):
        E+=[M[i]]
    return E

#INPUT : M    'a matrix
#        ind  'b list : tuple list 
#OUTPUT: List 'b list : support of column matrix M in tuple list
def Supp (M,ind):
    l=M.nrows()
    List=[]
    for i in range(0,l):
        if M[i]!=0:
            List+=[ind[i]]
    return List

#INPUT : M        'a matrix
#        A        'b list : tuple set
#        B        'b list : tuple set
#        Lignes   'b list : indices lines of M with tuples
#        Colonnes 'b list : indices columns of M with tuples
#OUTPUT: cMe      'a matrix : extracted matrix of M wrt A and B
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
            jc+=1
        ic+=1
    return matrixe(cMe)

#INPUT : M        'a matrix
#        A        'b list : tuple set
#        Sigma    'b list : tuple set
#        n        integer
#        Lignes   'b list        : indices lines of M with tuples
#        Colonnes 'b list        : indices columns of M with tuples
#OUTPUT: L        ('a list) list : n first independant lines of M
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

#3°)Tensorial product
#~~~~~~~~~~~~~~~~~~~~

#INPUT : M1  'a matrix : length m1*n1
#        m1  integer
#        n1  integer
#        M1  'a matrix : length m2*n2
#        m2  integer
#        n2  integer
#OUTPUT: res 'a matrix : length m1m2*n1n2 : tensorial product of M1 and M2
def Tproduit(M1,m1,n1,M2,m2,n2):
    return (MatrixSpace(QQ,m1*m2,n1*n2))((numpy.kron(M1,M2)).tolist())

#4°)Cartesian product
#~~~~~~~~~~~~~~~~~~~~

#INPUT : A 'a list 
#        B 'b list
#OUTPUT: E ('a * 'b) list : Cartesian product of A and B
def ProduitCart(A,B):
    E=[]
    i=0
    for ela in A:
        for elb in B:
            E+=[ela+elb]    
    return E

#5°)Main algorithm
#~~~~~~~~~~~~~~~~~

def ExtUpdate(Sigma,A,M1,M):
    #A priori le fait de numéroter par {0,1,2}, plutot que {-1,0,1} ne change rien 
    #et permet d'utiliser les memes fonctions que pour A .
    extSigma=ProduitCart([[0],[1],[2]],Sigma) 
    extA=ProduitCart([[0],[1],[2]],A)
    extM=Tproduit(M1,3,3,M,M.nrows(),M.ncols())
    return (extSigma,extA,extM)

def GrosCalcul(extA,P,dP,p,Q,deg,T,l,i,m):
    TaQi=Zero(len(extA))
    print("la longueur de extA est : %d" %len(extA))
    global PAR_SR
    if PAR_SR and a==3:
        ListArg=[]
        for Uplet in extA:
            R=dP
            r=p-1
            for j in range(i-1,m): 
                R*=(Q[j])**(Uplet[j-i+1])
                r+=deg[j]*(Uplet[j-i+1])
            ListArg+=[(l,T,R,r,P,p,Uplet)]
        Output=list(ParCalcul(ListArg))
        for j in range(len(extA)):
            Uplet=Output[j][0][0][6]
            TaQi[ConvertUplet(Uplet,extA)]=Output[j][1]
    else:
        for Uplet in extA :
            R=dP
            r=p-1
            for j in range(i-1,m): 
                R=R*(Q[j])**(Uplet[j-i+1])
                r=r+deg[j]*(Uplet[j-i+1])
            R,r=IntRem(l,T,R,r,P,p)
            TaQi[ConvertUplet(Uplet,extA)]=PmVPol(l,T,P,p,R,r)
    return TaQi

#INPUT : l     integer 
#        T     (int * Q[X1,...,Xl] * int) list : triangular system
#        P     Q[X1,...,Xl-1][Xl]
#        p     integer                 : degree of P
#        Q     Q[X1,...,Xl-1][Xl] list : list of derivatives of a polynomial QQ
#        deg   int list                : degrees of Q
#OUTPUT: Sigma (int list) list         : sign realization of QQ on roots of P
#        nb    (int list) list         : multiplicity order of roots
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
        extSigma,extA,extM=ExtUpdate(Sigma,A,M1,M)
        Mcoli=MatrixSpace(QQ,len(extA),1)
        TaQi=GrosCalcul(extA,P,dP,p,Q,deg,T,l,i,m)
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
