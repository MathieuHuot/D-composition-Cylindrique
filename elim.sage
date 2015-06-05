#(***************************************************************************************)
#(*                                                                                     *)
#(*                                                                                     *)
#(*                              HUOT Mathieu                                           *)
#(*                              GARNIER Remy                                           *)
#(*                    Licence 3 : stage de Mathématiques                               *)
#(*                           Version Q[X1,...Xn]                                       *)
#(*                           Phase d'élimination                                       *)
#(*                                                                                     *)
#(*                                                                                     *)
#(***************************************************************************************)

attach("sub_resultants.sage")
attach("int_rem.sage")
attach("pmv.sage")
attach("root_coding.sage")
attach("sign_degree.sage")
attach("sign_realization.sage")

#INPUT : P Q[X1,...,Xi]
#        i integer
#OUTPUT: Booléen : True ssi P n'est pas dans Q
def PasDansR(P,i):
    if i==0 or P==0:
        return false
    else:
        p=P.degree()
        if p > 0 :
            return true
        else:
            lcofP=P[p]
            return PasDansR(lcofP,i-1)
#COMPLEXITY : O(i)

#INPUT : l   integer
#        P   Q[X1,...,Xl-1][Xl]
#OUTPUT: res Q[X1,...,Xl-1][Xl] list : polynomes de la troncature de P
def Tru(l,P):
    Q=0*TdV[l-1]+P
    p=Q.degree()
    res=[]
    for r in range(p,-1,-1):
        if Q[r]!=0:
            res=res+[Q]
        if not PasDansR(P[r],l-1): #ie a_r est dans R*
            break
        Q=Q-Q[r]*TdV[l-1]**r
    return res
#COMPLEXITY : O(deg(P)*(l+deg(P)))

#INOUT : Q Q[X_1,...X_n-1][X_n] list
#OUTPUT: P Q[X_1,...X_n-1][X_n] list list
#        Construit les ensembles (P_i)_{i<=n} de la phase d'élimination
def Elim(Q):
    n=len(Q)
    P=Zero(n)  #On va diviser chaque polyome par son contenu et vérifier qu'il est dans
    P[n-1]=[Primitif(n,Pol) for Pol in Q[n-1]] #le bon anneau

    for i in range(n-1,-1,-1): #Construction inductive des P_i
        lon=len(P[i])
        B=TdA[i+1] #On se place dans l'anneau Q[X1][X2]...[X_i+1]
        for j in range(lon): #On rend P "séparable"
            P1=P[i][j]
            P2=diff(P1,TdV[i])
            P1=B(A(P1)//gcd(A(P1),A(P2)))
        for j in range(lon): #Si P_i divise P_j on peut simplifier en gardant
            for k in range(j):       #P_j/P_i et P_i
                P1=P[i][j]
                P2=P[i][k]
                P3=B(gcd(A(P1),A(P2)))
                p1=P1.degree()
                p2=P2.degree()
                p3=P3.degree()
                if p3>=p1 and IntRem2(i+1,P3,p3,P1,p1)==0:
                    P[i][k]=B(A(P[i][k])//A(P3))
                elif p3>=p2 and IntRem2(i+1,P3,p3,P2,p2)==0:
                    P[i][j]=B(A(P[i][j])//A(P3))
        NewA=[] #On enlève les polynomes devenus constants ou identiques
        for j in range(lon):
            on_veut_pas=false
            if not PasDansR(P[i][j]+0*TdV[i],i+1): #i.e le polynome est constant
                on_veut_pas=true
            for k in range(j+1,lon): #On regarde si on retrouve notre polynome plus loin
                if P[i][j]==P[i][k]:
                    on_veut_pas=true
                if on_veut_pas:
                    break
            if not on_veut_pas: #Alors on garde le polynome pour la suite
                NewA=NewA+[B(P[i][j])]
        P[i]=NewA #On met à jour après notre simplification

        if i==0:
            break
        P[i-1]=[Primitif(i,Pol) for Pol in Q[i-1]] # On initie P_i à Q_i
            
        for Pol in P[i]: # On calcule Elim_(X_i+i)(P_i)
            TruPol=Tru(i+1,Pol)
            for R in TruPol:
                r=R.degree()
                if PasDansR(R[r],i-1): #Si le coefficient dominant n'est pas dans QQ)
                    R2=A(R[r])       #On rend le polynome (qui est le LCoef) primitif
                    d=R2.content()
                    P[i-1]=P[i-1]+[TdA[i](R[r]*(1/d))] #Et on l'ajoute à P_i
                if r>1:
                    Sres=SubResultants(i+1,R,r,diff(R,TdV[i]),r-1)
                    for j in range(len(Sres)): #On rajoute les Sres_j non nuls
                        if PasDansR(0*TdV[i-1]+Sres[j],i):
                            Sres2=A(Sres[j])  #On rend le polynome primitif
                            d=Sres2.content()
                            P[i-1]=P[i-1]+[TdA[i](Sres[j]*(1/d))]
                for Pol2 in P[i]:
                    TruPol2=Tru(i+1,Pol2)
                    for T in TruPol2:
                        t=T.degree()
                        if t<=r:
                            if t==r: #On fait appel à Intrem avant SubResultants
                                T2=IntRem2(i+1,R,r,T,t)
                                if T2==0:
                                    Sres=[]
                                else:
                                    T2=B(T2*1/A(T2).content())
                                    t2=T2.degree()
                                    Sres=SubResultants(i+1,R,r,T2,t2)
                            else:
                                Sres=SubResultants(i+1,R,r,T,t)
                            for j in range(len(Sres)):
                                #On rajoute ceux qui ont un degré >0
                                if PasDansR(0*TdV[i-1]+Sres[j],i):
                                    d=A(Sres[j]).content() #On rend le polynome primitif
                                    P[i-1]=P[i-1]+[TdA[i](Sres[j]*(1/d))]
    return P
#COMPLEXITY :
