#(***************************************************************************************)
#(*                                                                                     *)
#(*                                                                                     *)
#(*                              GARNIER Remy                                           *)
#(*                              HUOT Mathieu                                           *)
#(*                    Licence 3 : stage de Mathématiques                               *)
#(*                           Version Q[X1,...Xn]                                       *)
#(*                             Line Partition                                          *)
#(*                                                                                     *)
#(*                                                                                     *)
#(***************************************************************************************)

#INPUT : a   int list
#        b   int list
#OUTPUT: res {-1,0,1} : 1 if a<b, 0 if a=b, -1 if b>a wrt Thom-encoding order
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
#NOTE  : Transforms a list of objects into a lists of singletons containing 
#        those objects and add the number or the root of P

def Singleton(L,P):
    res=[]
    l=len(L)
    for i in range(0,l):
        res+=[[[i+1,L[i],P]]] #we code a triplet with a list
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
#OUTPUT: res (int * int list * Q[X1,...,Xl]) list list optimized fusion sort of L1 and L2 
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
#OUTPUT: L  (int * int list * Q[X1,...,Xl]) list list : fusion of the lists of SL with TriP, removing
#           repetitions, i.e sort roots of the different polynomials with the help of their Thom-Encoding
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
#OUTPUT: Res Q[X1,...,Xl-1][Xl] P without its nul coefficients when taken in \alpha_1,...,\alpha_l-1
#        p   integer
def Normalize(l,T,P):
    p=Degree(l,T,P)
    Res=0
    for j in range(0,p+1):
        Res+=P[j]*TdV[l-1]**j #Main variable of P is X_l
    Res=Primitif(l,Res)    #We make it primitive
    return Res,p
#COMPLEXITY: O(2EXP)

#INPUT : PP2    Q[X1,...,Xl] list 
#        l      integer 
#        T      (int * Q[X1,...,Xl] * int)  list : triangular system
#OUTPUT: SL     (int * int list * Q[X1,...,Xl]) list list 
#        Normed (Q[X1,...,Xl] * int) list : normalized polynomials of PP2
#NOTE  : SL represents the real line partition by being a list of sample root codings 
#        of every intervalle of the partition
def LinePartition(PP2,l,T):
    lon2=len(PP2) #PP2 contains polynomials with l variables
    
    Normed2=Zero(lon2)
    if NORM_PAR:
        Par=[] #list to be parallelized
        for i in range(lon2):
            Par=Par+[(l,T,PP2[i],i)]
        Output=list(NormalizePar(Par))
        for i in range(lon2):
            j=Output[i][1][1]
            Normed2[j]=Output[i][1][0]
    else:
        Normed2=[Normalize(l,T,PP2[i]) for i in range(lon2)]
        
    Normed=[] #normalized polynomials with roots will go there
    RootCodePi=[] #and their rootcoding there
    if ROOT_PAR:

        List=[]
        for i in range(lon2):
            Pi,pi=Normed2[i]
            if pi>0:
                List=List+ [(l,T,Pi,pi,Pi,pi,i)]
        lon=len(List)
        Output=list(RootPar(List)) # Parallelisation of rootcodings
        for i in range(lon):
            k=0
            for j in range(lon):
                if Output[j][1][1]==i:
                    k=j
                    break
            Root=Output[k][1][0]
            if Root!=[]: #we only keep polynomials that do have roots
                Normed=Normed+[Normed2[i]]
                RootCodePi=RootCodePi+[Root]
    else:
        for i in range(lon2):
            Pi,pi=Normed2[i]
            if pi>0: #We dont want constant polynomials
                Root=RootCoding(l,T,Pi,pi,Pi,pi)           
                if Root!=[]: #we only keep polynomials that do have roots
                    Normed=Normed+[Normed2[i]]
                    RootCodePi=RootCodePi+[Root]
            
    lon=len(Normed)#We put all the root-codings of normalized polynomials in a mtatrix
    SLL=[[0 for j in range(lon)] for i in range(lon)]
    
    if ROOT_PAR2:
        ListArg=[] #List of the arguments for parallelization
        for i in range(lon):
            for j in range(lon):
                if i==j: #We already calculated and stocked the Root-Coding of Pi on its roots
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
                if i==j: #We already calculated and stocked the Root-Coding of Pi on its roots
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
            #We need SL[j] and j to know whether the root is a root of Q
        for j in range(i-1,-1,-1):
            SL[j]=EnlargeWithDroite(SL[j],SLL[j][i],Normed[i],SL[i],i)
    for i in range(lon):
        for j in range(len(SL[i])):
            SL[i][j]=[i+1]+SL[i][j] 
        #We write in the head of every coding of root the number of the polynomal that generates it
    SL=OrderedMerge(SL)
    return SL,Normed
#COMPLEXITY : O(2EXP)
