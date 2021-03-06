#(***************************************************************************************)
#(*                                                                                     *)
#(*                                                                                     *)
#(*                              GARNIER Remy                                           *)
#(*                              HUOT Mathieu                                           *)
#(*                    Licence 3 : stage de Mathématiques                               *)
#(*                           Version Q[X1,...Xn]                                       *)
#(*                             Int Reminder                                            *)
#(*                                                                                     *)
#(*                                                                                     *)
#(***************************************************************************************)

#INPUT : l   integer
#        T   (int * Q[X_1,...X_l] * int) list : triangular system
#        Q   Q[X_1,...X_l]
#        q   integer
#        P   Q[X_1,...X_l]
#        p   integer 
#OUTPUT: Res Q[X_1,...X_l] : polynomial proprotionnal to the rest in the euclidean divison of Q by P 
#        r   integer       : degree of res in T
def IntRem(l,T,Q,q,P,p):
    if q<p:
        return Q,q
    R=[Q[i] for i in [0..q+1]]
    for i in range(q-p,-1,-1):
        for j in range(0,p):
            R[i+j]=P[p]*R[i+j]-P[j]*R[i+p]
        for j in range(i):
            R[j]=P[p]*R[j]
    if (q-p)%2==0:
        for j in range(p):
            R[j]=R[j]*P[p]
    for i in range(p,q+1):
        R[i]=0
    Aux=InitPoly(l,R,q)
    r=Degree(l,T,Aux)
    Res=0
    for i in range(r+1):
        Res+=Aux[i]*TdV[l-1]**i
    #Res=Primitif(l,Res)#We make it primitive
    return Res,r
#COMPLEXITY : O(2EXP)
