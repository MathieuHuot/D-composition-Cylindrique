# -*-coding:Latin-1 -*

#(***************************************************************************************)
#(*                                                                                     *)
#(*                                                                                     *)
#(*                               GARNIER Remy                                          *)
#(*                               HUOT MAthieu                                          *)
#(*                    Licence 3 : stage de Mathématiques                               *)
#(*                            Version Q[X1,...Xn]                                      *)
#(*                            Fonctions générales                                      *)
#(*                                                                                     *)
#(*                                                                                     *)
#(***************************************************************************************)

#General utilitary functions
#===========================

#INPUT : n integer
#OUTPUT: s {-1,1} : signature of the permutation of n raws
def Epsilon(n):
    return (-1)**(n*(n-1)/2)
#COMPLEXITY : O(1)

#INPUT : (a,b) 'a * 'b
#OUTPUT: a     first element
def Fst(a,b):
    return a
#COMPLEXITY : O(1)

#INPUT : (a,b) 'a * 'b
#OUTPUT: b     second element
def Snd(a,b):
    return b
#COMPLEXITY : O(1)

#INPUT : n   integer
#OUTPUT: res int list : list of length n fulled with 0
def Zero(n):
    return [0 for i in range(0,n)]
#COMPLEXITY : O(n)

#INPUT : r rationnal
#OUTPUT: s {-1,0,1}   : sign of r  
def Sign_1(r):
    if r>0:
        return 1
    else: 
        if r<0:
            return -1
        else:
            return 0
#COMPLEXITY : O(1)

#INPUT : x   integer
#OUTPUT: x-1 integer
Decr= lambda x: x-1
#COMPLEXITY : O(1)

#INPUT : liste int list
#OUTPUT: liste int list : each element has been decreased by 1
def DecrL (liste):
    return map(Decr,liste)
#COMPLEXITY : O(len(liste))

#Polynomials rings
#=================

#Rings we currently use
A0=QQ
A1.<X1>=QQ[]
A2.<X2>=A1[]
A3.<X3>=A2[]
A4.<X4>=A3[]
A5.<X5>=A4[]
A6.<X6>=A5[]
A7.<X7>=A6[]
A8.<X8>=A7[]
A9.<X9>=A8[]
A10.<X10>=A9[]
A11.<X11>=A10[]
A12.<X12>=A11[]
A13.<X13>=A12[]
A14.<X14>=A13[]
A15.<X15>=A14[]
A16.<X16>=A15[]
A17.<X17>=A16[]
A18.<X18>=A17[]
A19.<X19>=A18[]
A20.<X20>=A19[]

#Variables Table
TdV=[X1,X2,X3,X4,X5,X6,X7,X8,X9,X10,X11,X12,X13,X14,X15,X16,X17,X18,X19,X20]
#Rings Table
TdA=[A0,A1,A2,A3,A4,A5,A6,A7,A8,A9,A10,A11,A12,A13,A14,A15,A16,A17,A18,A19,A20]
#Polynomial ring over QQ with n (here 7) variables
A=PolynomialRing(QQ,'X',7)

#Utilitary functions over  Q[X_1,...,X_l]:
#=========================================

#INPUT : l   integer
#        P   Q[X1,...,Xl]
#        p   integer : degree of P
#OUTPUT: Res Q[X1,...,Xl] list : list of the successive derivatives of P
#        Deg int list          : list of the degrees of the successive derivatives of P
def ListOfDerivate(l,P,p):
    Res=[P]
    deg=[p]
    Q=Primitif(l,P)
    for i in range(p-1,-1,-1):
        Q=diff(Q,TdV[l-1])
        Q=Primitif(l,Q)
        Res+=[Q]
        deg+=[i]
    return Res,deg
#COMPLEXITY : O(p*p)

#INPUT : l integer
#        R Q[X1,...,X(l-1)] list
#        q integer : degree of Q
#OUTPUT: Q Q[X_1,...,X_l]
def InitPoly (l,R,q):
    Q=0
    for i in range(0,q+1):
        Q+=R[i]*TdV[l-1]**i
    Q=Primitif(l,Q) #We make Q primitive
    return Q
#COMPLEXITY : O(q*q)

#INPUT : l   integer
#        P   Q[X1,...,Xl]
#        i   rationnal
#OUTPUT: Res Q[X1,...,Xl] : Res=P(X+i)
def Decale(l,P,i):
    p=P.degree()
    Res=P[p]
    for j in range(p-1,-1,-1): #Horner like
        Res=Res*(TdV[l-1]+i)+P[j] #Main variable of P is X_l
    Res=Primitif(l,Res) #Make it primitive
    return Res
#Complexity : O(deg(P)**2)

#INPUT : l integer
#        P Q[X1,...,Xl]
#OUTPUT: Q Q[X1,...,Xl] primitive version of P
def Primitif(l,P):
    if P==0:
        return 0
    else:
        Q=A(P)
        return TdA[l]((1/Q.content())*Q)
#COMPLEXITY : O(deg(P)*log(max(coef(P))))
