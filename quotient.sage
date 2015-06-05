# -*-coding:Latin-1 -*

#(***************************************************************************************)
#(*                                                                                     *)
#(*                                                                                     *)
#(*                              HUOT Mathieu                                           *)
#(*                              GARNIER Remy                                           *)
#(*                    Licence 3 : stage de Mathématiques                               *)
#(*                           Version Q[X1,...Xn]                                       *)
#(*                          Fonctions de Quotient                                      *)
#(*                                                                                     *)
#(*                                                                                     *)
#(***************************************************************************************)

attach("fonctions_generales.sage")

#INPUT : l    integer
#        P    Q[X1,...,Xl]
#        Q    Q[X1,...,Xl]
#OUTPUT: P//Q Q[X1,...,Xl] : division exacte de P par Q
#WARNING: ne fonctionne que si la divison exacte est possible selon notre ordre
#         sur les variables, ie : X_l > X_l-1 > ... > X1
def Quotient(l,P,Q):
    def Quo(n,P,Q):
        P=TdA[n](P)
        Q=TdA[n](Q)	
        if n==0:
            return QQ(P)/QQ(Q)
        else:
            p=P.degree()
            q=Q.degree()
            if P==0:
                return 0
            elif p<q:
                raise ValueError("P n'est pas divisible par Q")
            else:
                lcofR=Quo(n-1,P[p],Q[q])
                P=P-lcofR*TdV[n-1]**(p-q)*Q
                restR=Quo(n,P,Q)
                return (lcofR*TdV[n-1]**(p-q)+restR)
    try:
        return Quo(l,P,Q)
    except ValueError as err:
        print("P n'est pas divisible par Q", "P=%s" %P, "Q=%s" %Q)
#COMPLEXITY : O(2**l * 2**deg(P))        
