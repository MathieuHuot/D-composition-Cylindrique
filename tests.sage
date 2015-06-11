#(***************************************************************************************)
#(*                                                                                     *)
#(*                                                                                     *)
#(*                            HUOT Mathieu                                             *)
#(*                            GARNIER Remy                                             *)
#(*                  Licence 3 : stage de Mathématiques                                 *)
#(*                               TESTS                                                 *)
#(*                                                                                     *)
#(*                                                                                     *)
#(***************************************************************************************)

Tests=False
if Tests:
    #Sign
    T=[[2, X1^2 - 2, 2], [2, X2^2 - 7, 2], [1, 25*X3^3 - 64, 3]]
    P=P=X1*X3**2-X2
    print(Sign(3,T,P))
    T=[[2,X1**2-X1-1,2],[1,(2*X1-1)*X2**2-1,2]]
    P1=X2**2-1
    print(Sign(2,T,P1))
    P2=P1-X1**2
    print(Sign(2,T,P2))
    P3=X2**4-1/5
    print(Sign(2,T,P3))
    #Elim
    Q=[[X1, X1^2 - X1 - 1], [X2, (2*X1 - 1)*X2^2 - 1, X2 + X1^2 - 5]]
    print(Elim(Q))
    #Crée un élément d'un Thom-encoding fictif
    def CreeThom(n):
    s=[random.randint(-1,1) for j in range(0,n-1)]+[1]
    return s
    #int -> int list
    #Crée un Thom-encoding fictif pour tester le tri
    def CreeEncoding(n,p):
        s=Zero(n)
        for i in range(0,n):
            s[i]=[random.randint(-1,1) for j in range(0,p-1)]+[1]
        return s
    #int * int -> (int list) list
    #Lifting
    R=Elim(Q) 
    C=Lifting(R,Q)
    #Renvoie un polynome de degré inférieur ou égal à deg,
    #dont les coefficients sont entre mini et maxi
    import random
    import time

    def RandomPol(l,deg=4,mini=-10,maxi=10):
        if l==0:
         return (random.randint(mini,maxi))
        else:
            Q=0
            for j in range(0,deg+1):
                Q = Q+RandomPol(l-1,deg,mini,maxi)*TdV[l-1]**j
            return Q
    #int * int * int -> Q[X]

    def test(n):
        L=Zero(n)
        for i in range(n):
            print("Test %d" %(i+1))
            L=[[],[],[(X3**3)*RandomPol(2,3,-2,2)]]
            t0=time.clock()
            C=Elim(L)
            t1=time.clock()
            print("Durée de l'Elim= %s" %(t1-t0))
            print("Nombre de polynome= %d" %(len(C[0])+len(C[1])))
            t2=time.clock()
            Lift=Lifting(C,L)
            print("Durée du Lifting=%s" %(t2-t1))

