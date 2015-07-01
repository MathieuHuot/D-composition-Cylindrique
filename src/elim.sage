#(***************************************************************************************)
#(*                                                                                     *)
#(*                                                                                     *)
#(*                              GARNIER Remy                                           *)
#(*                              HUOT Mathieu                                           *)
#(*                    Licence 3 : stage de Mathématiques                               *)
#(*                           Version Q[X1,...Xn]                                       *)
#(*                            Elimination part                                         *)
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
#OUTPUT: Booléen : True iff P is not a rationnal number
def PasDansR(P,i):
    if i==0 or P==0:
        return false
    else:
        P=TdA[i](P)
        p=P.degree()
        if p > 0 :
            return true
        else:
            lcofP=P[p]
            return PasDansR(lcofP,i-1)
#COMPLEXITY : O(i)

#INPUT : l   integer
#        P   Q[X1,...,Xl-1][Xl]
#OUTPUT: res Q[X1,...,Xl-1][Xl] list : polynomes of the truncature of P
def Tru(l,P):
    P=TdA[l](P)
    p=P.degree()
    res=[]
    for r in range(p,-1,-1):
        if P[r]!=0:
            res=res+[P]
        if not PasDansR(P[r],l-1): #ie a_r is in R*
            break
        P-=P[r]*TdV[l-1]**r
    return res
#COMPLEXITY : O(deg(P)*(l+deg(P)))

#INPUT : lis Q[X1,...,Xi] list
#        i   integer           : number of variables of the polynomials of lis
#OUTPUT: lis Q[X1,...,Xi] list : every polynomial became separable 
def Separable(lis,i):
    lon=len(lis)
    for j in range(lon): 
        P1=lis[j]
        P2=diff(P1,TdV[i-1])
        lis[j]=Quotient(i,P1,gcd(A(P1),A(P2)))
    return lis
#COMPLEXITY : O(lon * 2**i * 2**(degmax(P_i)) )

#INPUT : lis  Q[X1,...,Xi] list
#        i    integer           : number of variables of the polynomials of lis
#OUTPUT: NewA Q[X1,...,Xi] list : the list that is simplified
def Simplify_1(lis,i):
    lon=len(lis)
    for j in range(lon): #If P_i divides P_j we may simplify by keeping
        for k in range(j):       #P_j//P_i and P_i
            P1=TdA[i](lis[j])
            P2=TdA[i](lis[k])
            P3=TdA[i](gcd(A(P1),A(P2)))
            p1=P1.degree()
            p2=P2.degree()
            p3=P3.degree()
            if IntRem2(i,P3,p3,P1,p1)==0:
                lis[k]=Quotient(i,P2,P3)
            elif IntRem2(i,P3,p3,P2,p2)==0:
                lis[j]=Quotient(i,P1,P3)
    NewA=[] #We remove polynomials that are the same or constant
    for j in range(lon):
        on_veut_pas=false
        if not PasDansR(lis[j],i): #i.e the polynomial is constant
            on_veut_pas=true
        for k in range(j+1,lon): #We look if we may find the polynomial further
            if lis[j]==lis[k] or lis[j]==-lis[k]:
                on_veut_pas=true
            if on_veut_pas:
                break
        if not on_veut_pas: #Then we keep the polynomial
            NewA=NewA+[lis[j]]
    return NewA
#COMPLEXITY : O(lon**2 * 2**i * degmax(P_i) )

#INPUT : P Q[X1,...,Xn] list list
#        i integer
#OUTPUT: P Q[X1,...,Xn] list list : P[i] is maybe simplified 
def Simplify_2(P,i):
    lis=P[i-1]
    NewA=[] #We are going to reduce the degree of polynomials by dividing p_j and p_(j+1)
    #by their gcd and by adding this gcd to the list if not constant
    lon=len(lis)
    for j in range(lon-1):
        P1=lis[j]
        P2=lis[j+1]
        P3=gcd(A(P1),A(P2))
        if PasDansR(P3,i): #Non trivial gcd
            P3=Primitif(i,P3) #We make it primitive
            lis[j]=Quotient(i,P1,P3) #We simplify P_j and P_j+1
            lis[j+1]=Quotient(i,P2,P3)
            if P3.degree()>0: #If it's of non nul degree in X_i+1 we add it directly
                NewA+=[P3]
            else: #Else we add it in the good level by watching for which variable its degree is >0
                k=i  
                while (TdA[k](P3)).degree()==0:
                    k-=1
                P[k-1]+=[P3]
        if lis[j].degree()>0: #If p_i is always of degree >=1, we keep it
            NewA+=[lis[j]]
    if lis[lon-1].degree()>0: #We add the last polynomial to the list
        NewA+=[lis[lon-1]]
    NewB=[] #We remove again polynomials that are now constant
    lon2=len(NewA)
    for j in range(lon2):
        on_veut_pas=false
        if not PasDansR(NewA[j],i): #i.e the polynomial is constant
            on_veut_pas=true
        for k in range(j+1,lon2): #We catch whether we find the polynomial further
            if NewA[j]==NewA[k] or NewA[j]==-NewA[k]:
                on_veut_pas=true
            if on_veut_pas:
                break
        if not on_veut_pas: #Then we keep the polynomial
            NewB=NewB+[NewA[j]]
        P[i-1]=NewB #We update after our simplifications
    return P 
#COMPLEXITY : O(lon * 2**i * 2**(degmax(P_i)) )

#INOUT : Q Q[X_1,...X_n-1][X_n] list
#OUTPUT: P Q[X_1,...X_n-1][X_n] list list
#        Build the sets (P_i)_{i<=n} of the elim phase
def Elim(Q):
    n=len(Q)
    P=[[] for i in range(n)]  #We will divide every polynomial by its content and 
    P[n-1]=[Primitif(n,Pol) for Pol in Q[n-1]] #check it's in the good ring

    for i in range(n-1,-1,-1): #inductive construction of P_i
        P[i]=Separable(P[i],i+1) #We make every polynomial separable
        P[i]=Simplify_1(P[i],i+1)
        P=Simplify_2(P,i+1)
        if i==0:
            break
        P[i-1]=[Primitif(i,Pol) for Pol in Q[i-1]] #We initialize P_i to Q_i
            
        for Pol in P[i]: #We calculate Elim_(X_i+i)(P_i)
            TruPol=Tru(i+1,Pol)
            for R in TruPol:
                R=TdA[i+1](R)
                r=R.degree()
                if PasDansR(R[r],i): #If leading coef is not in QQ we add it
                    P[i-1]=P[i-1]+[Primitif(i,R[r])]
                if r>1:
                    Sres=SubResultants(i+1,R,r,diff(R,TdV[i]),r-1)
                    for j in range(len(Sres)): #We add non nul Sres_j
                        if PasDansR(Sres[j],i):
                            P[i-1]=P[i-1]+[Primitif(i,Sres[j])]
                for Pol2 in P[i]:
                    TruPol2=Tru(i+1,Pol2)
                    for T in TruPol2:
                        T=TdA[i+1](T)
                        t=T.degree()
                        if t<=r:
                            if t==r: #We make a call to IntRem before SubResultant
                                T2=IntRem2(i+1,R,r,T,t)
                                if T2==0:
                                    Sres=[]
                                else:
                                    T2=Primitif(i+1,T2)
                                    t2=T2.degree()
                                    Sres=SubResultants(i+1,R,r,T2,t2)
                            else:
                                Sres=SubResultants(i+1,R,r,T,t)
                            for j in range(len(Sres)):
                                #We add those of degree >0
                                if PasDansR(Sres[j],i):
                                    P[i-1]=P[i-1]+[Primitif(i,Sres[j])]
    return P
#COMPLEXITY : O(2EXP)
