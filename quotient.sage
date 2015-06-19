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
    F=A(Q)
    Q=0
    while P.lm()!=0:
        D=P.lt()//F.lt()
        Q+=D
        P-=D*F
    return TdA[l](Q)
#COMPLEXITY : O(p**(l+1))        
