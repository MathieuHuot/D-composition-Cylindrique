#(***************************************************************************************)
#(*                                                                                     *)
#(*                                                                                     *)
#(*                              GARNIER Remy                                           *)
#(*                              HUOT Mathieu                                           *)
#(*                    Licence 3 : stage de Mathématiques                               *)
#(*                           Version Q[X1,...Xn]                                       *)
#(*                          Exact divison function                                     *)
#(*                                                                                     *)
#(*                                                                                     *)
#(***************************************************************************************)

attach("fonctions_generales.sage")

#INPUT : l    integer
#        P    Q[X1,...,Xl]
#        Q    Q[X1,...,Xl]
#OUTPUT: P//Q Q[X1,...,Xl] : exact divison of P by Q
#WARNING: works iff exact division is possible wrt usual order on variable
#         ie : X_l > X_l-1 > ... > X1
def Quotient(l,P,Q):
    return TdA[l](A(P)//A(Q))
#COMPLEXITY : O(p**(l+1))        
