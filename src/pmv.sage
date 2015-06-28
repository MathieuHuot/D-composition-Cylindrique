#(***************************************************************************************)
#(*                                                                                     *)
#(*                                                                                     *)
#(*                              GARNIER Remy                                           *)
#(*                              HUOT Mathieu                                           *)
#(*                    Licence 3 : stage de Mathématiques                               *)
#(*                           Version Q[X1,...Xn]                                       *)
#(*                      Permanence Minus Variations                                    *)
#(*                                                                                     *)
#(*                                                                                     *)
#(***************************************************************************************)

#INPUT : s   int list
#OUTPUT: res integer : Permanence Minus Variation of s
def PmV(s):
    res=0
    i=0
    count=1
    a=len(s)
    while a-i>1:
        while (i+count<a and s[i+count]==0):
            count+=1
	if i+count==a:
	   return res
        if ((count%2)!=0):
            res+=Epsilon(count)*Sign_1(s[i]*s[i+count])
        i+=count
        count=1
    return res
#COMPLEXITY : O(len(s))

#INPUT : l   integer
#        T   (int * Q[X_1,...X_l] * int) list : triangular system of level l
#        P   Q[X_1,...X_l]
#        p   integer
#        Q   Q[X_1,...X_l]
#        q   integer
#OUTPUT: res int      : PMV of the sign valuation of the sub-resultant of P and Q in T
def PmVPol(l,T,P,p,Q,q):
    s=Zero(p+1)
    if q==-1:
        return 0
    sRes=SubResultants(l,P,p,Q,q)
    for j in range(0,p+1):
        s[j]=Sign(l-1,T,sRes[j]) #l-1 because Sres[j] are in Q[X_1,...,X_l-1]
    return PmV(s)
#COMPLEXITY : O(2EXP) 
