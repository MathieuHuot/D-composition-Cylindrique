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
def Quotient(l,P,Q,numero=50):
    P=A(P)
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
            else:
                lcofR=Quo(n-1,P[p],Q[q])
                P=P-lcofR*TdV[n-1]**(p-q)*Q
                restR=Quo(n,P,Q)
                return (lcofR*TdV[n-1]**(p-q)+restR)
    if len(P.monomials())>numero:
        return Quo(l,P,Q)
    else:
        F=A(Q)
        Q=0
        while P.lm()!=0:
            D=P.lt()//F.lt()
            Q+=D
            P-=D*F
        return TdA[l](Q)
#COMPLEXITY : O(2**l * 2**deg(P))        
