#(***************************************************************************************)
#(*                                                                                     *)
#(*                                                                                     *)
#(*                              GARNIER Remy                                           *)
#(*                              HUOT Mathieu                                           *)
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
#NOTE : More basic version which is an EnlargeWithGauche with no root searching
def EnlargeWithCompleting(lis,code,Q):
    n=len(lis)
    res=[]
    for i in range(n):
        res=res+[[[-1,code[i],Q]]+lis[i]]            
    return res
#COMPLEXITY : O(len(lis))

#INPUT : shortL (int * int list * Q[X_1,...,X_l]) list : a list of coding of roots
#        v      int list : Thom Encoding of P for a root \alpha_i
#        oldv   int list : Thom Encoding of oldP for a root \alpha_j < \alpha_i
#        P      Q[X_1,...,X_l]
#        oldP   Q[X_1,...,X_l]
#OUTPUT: F      int * (int * int list * Q[X_1,...,X_l]) list 
#               : root in shortL (vp,P) such that vp<v, voldP>oldv for the Thom Encoding order
def Find(shortL,v,oldv,P,oldP):
    n=len(shortL)
    r1=shortL[0]
    m=len(r1)
    indiceP=-1
    indiceOP=-1
    for i in range(m): #We find i and j indices such that P=Pi, oldP=Pj
        if P==r1[i][2]:
            indiceP=i
            break
    for i in range(m):
        if oldP==r1[i][2]:
            indiceOP=i
            break
    for i in range(n): #We look for the wanted root F among roots L[i] in L
        F=shortL[i]
        vP=F[indiceP][1]
        voldP=F[indiceOP][1]
        if PetitThom(vP,v) and PetitThom(oldv,voldP):
            Ind=len(F) #Same technique as before: 1st element is the indice of F where
            return [Ind]+F # "r is defined"
#COMPLEXITY : O(len(shortL) + len(shortL[0])) 

#INPUT : l  integer
#        T  (int * Q[X1,...,Xl] * int) list 
#        PP (Q[X1,...,Xl] * int) list 
#        i  integer
#        j  integer
#OUTPUT: shortL (int * int list * Q[X1,...,Xl]) list list : R-rootcoding of roots of (PiPj)'
def PreCalculCompleting(l,T,PP,i,j):
    m=len(PP)
    P,p=PP[i]
    Q,q=PP[j]
    pol=diff(P*Q,TdV[l-1])
    pol=Primitif(l,pol) #Make it primitive
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
#OUTPUT: shortL (int * int list * Q[X1,...,Xl]) list list : R-rootcoding of roots of Pi'
def PreCalculCompleting2(l,T,PP,i):
    m=len(PP)
    P,p=PP[i] #polynomial Pi and its degree
    P=diff(P,TdV[l-1]) #Pi'
    p-=1 #and its degree
    P=Primitif(l,P) #We make it primitive
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
#OUTPUT: shortL (int * int list * Q[X1,...,Xl]) list list : R-rootcoding of roots of (PoldP)'
def CalculCompleting(l,T,PP,oldP,P):
    m=len(PP)
    pol=diff(P*oldP,TdV[l-1])
    pol=Primitif(l,pol) #We make it primitive
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
#NOTE : Complete the real line of level l which is split by LinePartition by adding a 
#root to represent each cell which is not a singleton(singletons being roots of polynomials
#obtained during the elim phase

def Completing(l,T,L,PP):
    n=len(L) 
    m=len(PP)
#There are n roots so 2n+1 representants of cells which are to be stocked in newL
    newL=Zero(2*n+1) 
#We will pre-calcul root-codings of (PiPj)' where L[i]=root of Pi or Pj and L[i+1]=root of
#Pi or Pj for at least two different i. indices stocked in Mtemp and resuts in M
    M=[["Vide" for i in range(m)] for j in range(m)] #to stock Root-codings
    Mtemp=[[0 for i in range(m)] for j in range(m)] #to count if a pre-calculus is worth
    
    for i in range(n-1):
        e=L[i]
        f=L[i+1]
        Ind1=e[0]-1 #P_Ind1(e)=0
        Ind2=f[0]-1 #P_Ind2(f)=0
        Mtemp[Ind1][Ind2]+=1
    if PAR_PCC:
        MPar=[] #in order to parallelize pre-calculus
        for i in range(m): 
            for j in range(i):
                Mtemp[i][j]+=Mtemp[j][i]
                if Mtemp[i][j]>1: #worth to make pre-calculus
                    MPar+=[(l,T,PP,i,j)] 
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
            if Mtemp[i][i]>0: #ie we want to avoid the calculus of (PiPi)'
                MPar+=[(l,T,PP,i)]
        ResPar=list(ParPreCalculCompleting2(MPar))
        for k in range(len(ResPar)):
            i=ResPar[k][0][0][3]
            M[i][i]=ResPar[k][1]
    else:
        for i in range(m):
            if Mtemp[i][i]>0:
                M[i][i]=PreCalculCompleting2(l,T,PP,i)
    
    for i in range(n): #We recopy what does not change : cells which are singletons
        newL[2*i+1]=L[i] #containing roots of the polynomials of PP

    E=L[0] #Case of the cell - infinity
    r,v,P=E[E[0]] #E[0] gives the indice of the polynomial from which E codes a root
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
                    
    newL[0]=[len(newL[0][0])]+newL[0][0] #We are interested in the 1st root
    oldv=v
    oldP=P

    for i in range(1,n): #Generla case of a cell between two roots
        E=L[i]
        r,v,P=E[E[0]]
        Ind1=max(L[i-1][0],L[i][0])-1
        Ind2=min(L[i-1][0],L[i][0])-1
        if Ind1==Ind2: #we look is the result has already been calculated
            newL[2*i]=M[Ind1][Ind2]
        elif Mtemp[Ind1][Ind2]>1: 
            newL[2*i]=M[Ind1][Ind2]
        else: #else we calculate it now
            newL[2*i]=CalculCompleting(l,T,PP,oldP,P)
        #we stock the good root i.e the one between L[i] and L[i+1]
        newL[2*i]=Find(newL[2*i],v,oldv,P,oldP)
        oldv=v
        oldP=P

    E=L[n-1] #Case of the cell + infinity
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
            
    #We are interested in the biggest root
    newL[2*n]=newL[2*n][len(newL[2*n])-1]
    newL[2*n]=[len(newL[2*n])]+newL[2*n] #We again add the "i" that shows the position
    return newL                          # where "r is defined"
#COMPLEXITY : O(2EXP)
