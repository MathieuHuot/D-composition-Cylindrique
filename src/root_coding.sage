#(***************************************************************************************)
#(*                                                                                     *)
#(*                                                                                     *)
#(*                              GARNIER Remy                                           *)
#(*                              HUOT Mathieu                                           *)
#(*                    Licence 3 : stage de Mathématiques                               *)
#(*                           Version Q[X1,...Xn]                                       *)
#(*                               Root Coding                                           *)
#(*                                                                                     *)
#(*                                                                                     *)
#(***************************************************************************************)

#INPUT : a int list
#        b int list
#OUTPUT: Booléen : Compares a to b wrt Thom-encoding order and returns True iff a<b
def PetitThom(a,b):
    if a==b:
        return False
    k=len(a)-1
    while a[k]==b[k]:
        k-=1
    if (a[k+1]==1 and a[k]<b[k]) or (a[k+1]==-1 and a[k]>b[k]):
        return True
    return False
#COMPLEXITY : WORSE : O(len(a)) AVERAGE : O(1)

#INPUT : L (int list) list
#OUTPUT: L (int list) list : sorted wrt Thom-Encoding order
def TriRapide(L):
    def trirap(L, g, d):
        pivot = L[(g+d)//2]
        i = g
        j = d
        while True:
            while PetitThom(L[i],pivot):
                i+=1
            while PetitThom(pivot,L[j]):
                j-=1
            if i>j:
                break
            if i<j:
                L[i],L[j]=L[j],L[i]
            i+=1
            j-=1
        if g<j:
            trirap(L,g,j)
        if i<d:
            trirap(L,i,d)
    g=0
    d=len(L)-1
    if d-g==0:
        return L
    trirap(L,g,d)
    return L
#COMPLEXITY : AVERAGE O(len(L)*log(len(L)))

#INPUT : sigma (int list) list
#        nb    (int list) list 
#OUTPUT: res   (int list) list : sigma expanded wrt multiplicity given by nb
def SigmaEtendu(sigma,nb):
    res=[]
    for i in range(0,len(nb)):
        res+=[sigma[i] for k in range(0,nb[i][0])]
    return res
#COMPLEXITY : O(len(sigma)*max(len(nb[i][0])))

#INPUT : l  integer
#        T  (int * Q[X_1,...,X_l-1] * int) list : triangular system
#        Pl Q[X_1,...,X_l]
#        pl integer        : degree of Pl
#        P  Q[X_1,...,X_l]
#        p  integer        : degree of P
#OUTPUT: s (int list) list : P-encoding of Pl roots in T
def RootCoding(l,T,Pl,pl,P,p):
    s=Zero(p+1)
    PP,deg=ListOfDerivate(l,P,p)
    a=SignRealization(l,T,Pl,pl,PP,deg)
    if a==[] or a==None:
        return []
    else:
        sigma,nb=a
        s=SigmaEtendu(sigma,nb)
        if s==[]:
            return []
        else:
            return TriRapide(s)
#COMPLEXITY: O(2EXP)
