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

#INPUT : lis Q[X1,...,Xi] list
#        i   integer nombre de variable des polynomes de lis
#OUTPUT: lis Q[X1,...,Xi] list : chaque polynome est devenu séparable 
def Separable(lis,i):
    lon=len(lis)
    for j in range(lon): 
        P1=lis[j]
        P2=diff(P1,TdV[i-1])
        lis[j]=TdA[i](A(P1)//gcd(A(P1),A(P2)))
    return lis
#COMPLEXITY : O(lon * 2**i * 2**(degmax(P_i)) )

#INPUT : lis  Q[X1,...,Xi] list
#        i    integer nombre de variable des polynomes de lis
#OUTPUT: NewA Q[X1,...,Xi] list : la liste de départ éventuellement restreinte 
def Simplify_1(lis,i):
    lon=len(lis)
    for j in range(lon): #Si P_i divise P_j on peut simplifier en gardant
        for k in range(j):       #P_j/P_i et P_i
            P1=lis[j]
            P2=lis[k]
            P3=TdA[i](gcd(A(P1),A(P2)))
            p1=P1.degree()
            p2=P2.degree()
            p3=P3.degree()
            if p3>=p1 and IntRem2(i,P3,p3,P1,p1)==0:
                lis[k]=TdA[i](A(lis[k])//A(P3))
            elif p3>=p2 and IntRem2(i,P3,p3,P2,p2)==0:
                lis[j]=TdA[i](A(lis[j])//A(P3))
    NewA=[] #On enlève les polynomes devenus constants ou identiques
    for j in range(lon):
        on_veut_pas=false
        if not PasDansR(lis[j],i): #i.e le polynome est constant
            on_veut_pas=true
        for k in range(j+1,lon): #On regarde si on retrouve notre polynome plus loin
            if lis[j]==lis[k] or lis[j]==-lis[k]:
                on_veut_pas=true
            if on_veut_pas:
                break
        if not on_veut_pas: #Alors on garde le polynome pour la suite
            NewA=NewA+[TdA[i](lis[j])]
    return NewA
#COMPLEXITY : O(lon**2 * 2**i * degmax(P_i) )

#INPUT : P Q[X1,...,Xn] list list
#        i integer
#OUTPUT: P Q[X1,...,Xn] list list : P[i] est éventuellement simplifié 
def Simplify_2(P,i):
    lis=P[i-1]
    NewA=[] #On va réduire le degré de polynômes en divisant P_j et P_(j+1) par leur pgcd et 
    #en ajoutant ce dernier à la liste s'il est non constant
    lon=len(lis)
    for j in range(lon-1):
        P1=lis[j]
        P2=lis[j+1]
        P3=TdA[i](gcd(A(P1),A(P2)))
        if PasDansR(P3,i): #On a un pgcd non trivial
            P3=Primitif(i,P3) #On le rend primitif
            lis[j]=TdA[i](A(P1)//A(P3)) #On simplifie P_j et P_j+1
            lis[j+1]=TdA[i](A(P2)//A(P3))
            if P3.degree()>0: #Si il est de degré non nul en X_i+1 on l'ajoute directement
                NewA+=[P3]
            else: #Sinon on l'ajoute dans le bon niveau en regardant pour quel variable son degré 
                k=i   #est strictement positif
                while (TdA[k](P3)).degree()==0:
                    k-=1
                P[k-1]+=[P3]
        if lis[j].degree()>0: #Si P_j est toujours de degré >=1, on le garde
            NewA+=[lis[j]]
    if lis[lon-1].degree()>0: #On rajoute le dernier polynome a la liste
        NewA+=[lis[lon-1]]
    NewB=[] #On re-enlève les polynomes devenus constants ou identiques
    lon2=len(NewA)
    for j in range(lon2):
        on_veut_pas=false
        if not PasDansR(NewA[j],i): #i.e le polynome est constant
            on_veut_pas=true
        for k in range(j+1,lon2): #On regarde si on retrouve notre polynome plus loin
            if NewA[j]==NewA[k] or NewA[j]==-NewA[k]:
                on_veut_pas=true
            if on_veut_pas:
                break
        if not on_veut_pas: #Alors on garde le polynome pour la suite
            NewB=NewB+[NewA[j]]
        P[i-1]=NewB #On met à jour après nos simplifications
    return P 
#COMPLEXITY : O(lon * 2**i * 2**(degmax(P_i)) )

#INOUT : Q Q[X_1,...X_n-1][X_n] list
#OUTPUT: P Q[X_1,...X_n-1][X_n] list list
#        Construit les ensembles (P_i)_{i<=n} de la phase d'élimination
def Elim(Q):
    n=len(Q)
    P=[[] for i in range(n)]  #On va diviser chaque polyome par son contenu et vérifier qu'il est
    P[n-1]=[Primitif(n,Pol) for Pol in Q[n-1]] #dans le bon anneau

    for i in range(n-1,-1,-1): #Construction inductive des P_i
        B=TdA[i+1] #On se place dans l'anneau Q[X1][X2]...[X_i+1]
        P[i]=Separable(P[i],i+1) #On rend chaque polynome à racine simple
        P[i]=Simplify_1(P[i],i+1)
        P=Simplify_2(P,i+1)
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
#COMPLEXITY : O(2EXP)
