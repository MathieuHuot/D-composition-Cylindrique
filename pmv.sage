#(***************************************************************************************)
#(*                                                                                     *)
#(*                                                                                     *)
#(*                              HUOT Mathieu                                           *)
#(*                              GARNIER Remy                                           *)
#(*                    Licence 3 : stage de Mathématiques                               *)
#(*                           Version Q[X1,...Xn]                                       *)
#(*                      Permanence minus variations                                    *)
#(*                                                                                     *)
#(*                                                                                     *)
#(***************************************************************************************)

#INPUT : s   int list
#OUTPUT: res integer : Permanence Minus Variation de s
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
#        T   (int * Q[X_1,...X_l] * int) list : système triangulaire de niveau l
#        P   Q[X_1,...X_l]
#        p   integer
#        Q   Q[X_1,...X_l]
#        q   integer
#OUTPUT: res int      : PMV de la valuation de signe du sous-résultant de P et Q dans T
def PmVPol(l,T,P,p,Q,q):
    s=Zero(p+1)
    if q==-1:
        return 0
    sRes=SubResultants(l,P,p,Q,q)
    for j in range(0,p+1):
        s[j]=Sign(l-1,T,sRes[j]) #l-1 car les Sres[j] sont dans Q[X_1,...,X_l-1]
    return PmV(s)
#COMPLEXITY : O(2EXP) 
