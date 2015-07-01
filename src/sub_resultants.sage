#(***************************************************************************************)
#(*                                                                                     *)
#(*                                                                                     *)
#(*                            GARNIER Remy                                             *)
#(*                            HUOT Mathieu                                             *)
#(*                  Licence 3 : stage de Mathématiques                                 *)
#(*                           Sub Resultants                                            *)
#(*                                                                                     *)
#(*                                                                                     *)
#(***************************************************************************************)

attach("quotient.sage")

#INPUT : l   integer
#        Q   Q[X_1,...X_l]
#        q   integer
#        P   Q[X_1,...X_l]
#        p   integer
#OUTPUT: Res Q[X_1,...X_l] : polynomial proportionnal to the rest in the division of Q by P
#Note  : Does not calculate the degree of Res in a triangular system, used in elim phase
def IntRem2(l,Q,q,P,p):
    if p>q:
        return Q
    R= [Q[i] for i in range(q+1)]
    for i in range(q-p,-1,-1):
        for j in range(0,p):
            R[i+j]= P[p]*R[i+j]-P[j]*R[i+p]
        for j in range(0,i):
            R[j]= P[p]*R[j]
    if (q-p)%2 == 0 :
        for j in range(0,p):
            R[j]= R[j]*P[p]
    for i in range(p,q+1):
        R[i]= 0
    Res = InitPoly(l,R,q) #Res est primitive
    return Res
#COMPLEXITY: O(q²-p²)

#INPUT : l integer
#        P Q[X_1,...X_l]
#        p integer
#        Q Q[X_1,...X_l]
#        q integer
#OUTPUT: s  Q[X_1,...X_l-1] list : sub-resultants of P and Q
def SubResultants(l,P,p,Q,q):
    s=Zero(p+1) # values in Q[X_1,...X_l-1]
    t=Zero(p+1) # values in Q[X_1,...X_l-1]
    SresP=Zero(p+1) # values in Q[X_1,...X_l-1][X_l]
    SresP[p]=P
    s[p]=1  
    t[p]=1  
    SresP[p-1]=Q
    t[p-1]=Q[q]
    if q == p-1 :
        s[p-1]=t[p-1]
    else:
        s[p-1]=0
    for li in range(q+1,p-1):
        s[li]=0
    SresP[q]=Epsilon(p-q)*(t[p-1])**(p-q-1)*Q
    t[q]=SresP[q][q]
    s[q]=t[q]
    i=p+1
    j=p
    while q!=-1:
        k=q
        t[j-1]=SresP[j-1][k]
        if k==j-1:
            s[j-1]=t[j-1]
        else :
            for li in range(k+1,j) :
                s[li]=0               #We know l>=1
            t[k]=Epsilon(j-k)*Quotient(l-1,(t[j-1])**(j-k),(s[j])**(j-k-1)) 
            s[k]=t[k]
            SresP[k]=Quotient(l,s[k]*SresP[j-1],t[j-1])
        Pol=IntRem2(l,t[j-1]*s[k]*SresP[i-1],j,SresP[j-1],k) #Pol is in Q[X_1,...X_l-1][X_l]
        if Pol==0:
            q=-1
        else:
            q=Pol.degree()
        div=s[j]*t[i-1]+0*TdV[l-1] #values in Q[X_1,...X_l-1][X_l]
        SresP[k-1]=-Quotient(l,Pol,div)
        i=j
        j=k
    for li in range(0,j-1):
        s[li]=0
    s[p]=P[p]
    #for i in range(p+1):
    #    s[i]=Primitif(l-1,s[i]) #We make it primitive
    return s
#COMPLEXITY : O(p*2^l*2^p)
