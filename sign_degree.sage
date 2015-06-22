#(***************************************************************************************)
#(*                                                                                     *)
#(*                                                                                     *)
#(*                            GARNIER Remy                                             *)
#(*                            HUOT Mathieu                                             *)
#(*                  Licence 3 : stage de Mathématiques                                 *)
#(*                           Sign & Degree                                             *)
#(*                                                                                     *)
#(*                                                                                     *)
#(***************************************************************************************)

#INPUT : l    integer
#        T    (int * Q[X_1,...X_l] * int) list : système triangulaire
#        P    Q[X_1,...X_l]
#OUTPUT: v[0] {-1,0,1} : sign de P en les racines donnés par le système triangulaire T
def Sign(l,T,P):
    if l==0:
        return Sign_1(QQ(P))
    p=Degree(l,T,P)
    if p==-1:
        return 0
    sigma=RootCoding(l,T,T[l-1][1],T[l-1][2],P,p)
    v=sigma[T[l-1][0]-1]  #car on numérote a partir de la 1ère racine et non 0ème
    return v[0]
#COMPLEXITY : O(2EXP)

#Renvoie le degré d'un polynome P dans un système triangulaire T de niveau l
#P est dans D[X],  D=Q[X_1,...X_l-1], X=X_l

#INPUT : l   integer
#        T   (int * Q[X_1,...X_l] * int) list : système triangulaire
#        P   Q[X_1,...X_l]
#OUTPUT: res integer : degré de P(Xl) dans le système triangulaire T
def Degree(l,T,P):
    if P in QQ:
        if P==0:
            return -1
        return 0
    i=0
    P=0*TdV[l-1] + P #s'assurer que P est vu dans Q[X_1,...X_l]
    p=P.degree() 
    s=Sign(l-1,T,P[p])
    while s==0:
        i+=1
        if i>p:
            return -1
        s=Sign(l-1,T,P[p-i])
    return (p-i)
#COMPLEXITY : O(2EXP)
