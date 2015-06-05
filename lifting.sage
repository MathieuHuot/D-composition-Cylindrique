#(***************************************************************************************)
#(*                                                                                     *)
#(*                                                                                     *)
#(*                              HUOT Mathieu                                           *)
#(*                              GARNIER Remy                                           *)
#(*                    Licence 3 : stage de Mathématiques                               *)
#(*                           Version Q[X1,...Xn]                                       *)
#(*                                Lifting                                              *)
#(*                                                                                     *)
#(*                                                                                     *)
#(***************************************************************************************)

test=False
debug=False

attach("elim.sage")
attach("line_partition.sage")
attach("completing.sage")
attach("print.sage")
if test:
    attach("tests.sage")

#INPUT : P   Q[X1,...,Xl]
#        rac integer * (int * int list * Q[X1,...,Xl]) list
#OUTPUT: i   integer : position de P dans rac si présent et 0 sinon
def RechP (P,rac):
    m=len(rac)
    notfound=true
    i=m-1
    while i > 0 and notfound:
        if rac[i][2]==P:
            notfound=false
        else: 
            i=i-1
    return i
#COMPLEXITY : O(len(rac))

#Crée l'arbre de la phase de remontée
#INPUT : PPtot  (Q[X1,...,Xn] list) list
#        PPlist (Q[X1,...,Xn] list) list
#OUTPUT: foret  (Q[X1,...,Xn] * int) tree : l'arbre de remontée contenant les polynomes et 
#                                           leur valuation de signe
def Lifting(PPtot,PPlist):
    k=len(PPtot)
    def Lift(l,T): #Construction récursive de chaque niveau
        L,PP=LinePartition(PPtot[l-1],l,T)
        m=len(PPlist[l-1])
        if L==[]: #Aucun polynome n'a de racine
            Tbis=T+[[1,TdV[l-1],1]] #X_l devient représentant de la ligne réelle
            Teval=[]
            for j in range(m):
                P=PPlist[l-1][j]
                p=P.degree()
                s=Sign(l-1,T,P[p])  #Le signe d'un polynome sans racine réelle
                Teval=Teval+[(P,s)] #est celui de son coefficient dominant
            if l<k:   #Si il reste un niveau à construire, on appelle récursivement
                arb=[Teval,Lift(l+1,Tbis)]
            else:     #Sinon c'est qu'on est arrivé à une feuille
                arb=[Teval,[]] 
            return [arb]
        else:
            foret=[]  #La ligne réelle est scindée par des racines de polynomes
            eval=[]   #On appelle donc completing pour avoir un représentant de
            L=Completing(l,T,L,PP) #de chaque cellule
            for i in range(len(L)):
                Teval=[]
                ind=L[i][0] #L'indice i d'un Pi tel que L[i] code une racine de Pi 
                P=L[i][ind][2]
                r=L[i][ind][0]
                Tbis=T+[(r,P,Degree(l,T,P))] #Qu'on ajoute au système triangulaire
                for j in range(m):
                    Pol=PPlist[l-1][j]
                    pos=RechP(Pol,L[i])
                    if pos==0:
                        deg=Pol.degree()
                        sig=Sign(l,Tbis,Pol)
                        Teval=Teval+[(Pol,sig)]
                    else:
                        Teval=Teval+[(L[i][pos][2],L[i][pos][1][0])]
                if l<k: #On appelle récursivement sur chaque noeud la construction du 
                    foret=foret+[[Teval,Lift(l+1,Tbis)]] #niveau suivant
                else: #Ou alors on est arrivé au plus bas niveau et on a des feuilles
                    foret=foret+[[Teval,[]]]
            return foret
    return [[],Lift(1,[])] #Construction à partir de la racine
#COMPLEXITY : O(2EXP)
