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
attach("parallelize.sage")
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

def Eval(L,T,l,PPlist,i):
    Teval=[]
    lon=len(PPlist[l-1])
    ind=L[i][0] #L'indice i d'un Pi tel que L[i] code une racine de Pi 
    P=L[i][ind][2]
    r=L[i][ind][0]
    Tbis=T+[(r,P,Degree(l,T,P))] #Qu'on ajoute au système triangulaire
    for j in range(lon):
        Pol=PPlist[l-1][j]
        pos=RechP(Pol,L[i])
        if pos>0:
            Teval=Teval+[(L[i][pos][2],L[i][pos][1][0])]
        else:
            fini=False #tant qu'on peut simplifier
            trouve=True #si on a trouvé une simplification : on refait une boucle
            sP=1 #Sign de Pol
            while (not fini) and trouve: 
                for m in range(1,len(L[i])):
                    trouve=False
                    Pol2=L[i][m][2]
                    Al=IntRem2(l,Pol,Pol.degree(),Pol2,Pol2.degree())
                    if Al==0:
                        trouve=True
                        if L[i][m][1][0]==0: #i.e Pol2(\alpha)=0
                            sP=0
                            fini=True
                            break
                        else:
                            Pol=Quotient(l,Pol,Pol2)
                            sP=sP*L[i][m][1][0]
            sP=sP*Sign(l,Tbis,Pol)
            Teval=Teval+[(Pol,sP)]
    return Teval,Tbis

def Eval2(L,T,l,PPlist,i,Tbis):
    Teval=[]
    lon=len(PPlist[l-1])
    for j in range(lon):
        Pol=PPlist[l-1][j]
        pos=RechP(Pol,L[i])
        if pos>0:
            Teval=Teval+[(L[i][pos][2],L[i][pos][1][0])]
        else:
            fini=False #tant qu'on peut simplifier
            trouve=True #si on a trouvé une simplification : on refait une boucle
            sP=1 #Sign de Pol
            while (not fini) and trouve: 
                for m in range(1,len(L[i])):
                    trouve=False
                    Pol2=L[i][m][2]
                    Al=IntRem2(l,Pol,Pol.degree(),Pol2,Pol2.degree())
                    if Al==0:
                        trouve=True
                        if L[i][m][1][0]==0: #i.e Pol2(\alpha)=0
                            sP=0
                            fini=True
                            break
                        else:
                            Pol=Quotient(l,Pol,Pol2)
                            sP=sP*L[i][m][1][0]
            sP=sP*Sign(l,Tbis,Pol)
            Teval=Teval+[(Pol,sP)]
    return Teval

#Crée l'arbre de la phase de remontée
#INPUT : PPtot  (Q[X1,...,Xn] list) list
#        PPlist (Q[X1,...,Xn] list) list
#OUTPUT: foret  (Q[X1,...,Xn] * int) tree : l'arbre de remontée contenant les polynomes et 
#                                           leur valuation de signe
def Lifting(PPtot,PPlist):
    k=len(PPtot)
    return [[],Lift(PPtot,PPlist,1,[],k)] #Construction à partir de la racine
#COMPLEXITY : O(2EXP)

def Lift(PPtot,PPlist,l,T,k): #Construction récursive de chaque niveau
    L,PP=LinePartition(PPtot[l-1],l,T)
    if L==[]: #Aucun polynome n'a de racine
        lon=len(PPlist[l-1])
        Tbis=T+[(1,TdV[l-1],1)] #X_l devient représentant de la ligne réelle
        Teval=[]
        for j in range(lon):
            P=PPlist[l-1][j]
            p=P.degree()
            s=Sign(l-1,T,P[p])  #Le signe d'un polynome sans racine réelle
            Teval=Teval+[(P,s)] #est celui de son coefficient dominant
        if l<k:   #Si il reste un niveau à construire, on appelle récursivement
            arb=[Teval,Lift(PPtot,PPlist,l+1,Tbis,k)]
        else:     #Sinon c'est qu'on est arrivé à une feuille
            arb=[Teval,[]] 
        return [arb]
    else:
        foret=[]  #La ligne réelle est scindée par des racines de polynomes
        #On appelle donc completing pour avoir un représentant de de chaque cellule
        L=Completing(l,T,L,PP) 
        if l==1:
            ListArg=[] #Début de la séquence de parellization
            for i in range(len(L)):
                ind=L[i][0] #L'indice i d'un Pi tel que L[i] code une racine de Pi 
                P=L[i][ind][2]
                r=L[i][ind][0]
                ListArg+=[(PPtot,PPlist,l+1,T+[(r,P,Degree(l,T,P))],k,i)]
            Output=list(LiftPar(ListArg))
            for i in range(len(L)):
                for j in range(i,len(L)):
                    if Output[i][1][1]>Output[j][1][1]:
                        Output[i],Output[j]=Output[j],Output[i]
            for i in range(len(L)):
                Teval=Eval2(L,T,l,PPlist,i,ListArg[i][3])
                foret+=[[Teval,Output[i][1][0]]]
        else:    
            for i in range(len(L)):
                Teval,Tbis=Eval(L,T,l,PPlist,i)
                if l<k: #On appelle récursivement sur chaque noeud la construction du 
                    foret+=[[Teval,Lift(PPtot,PPlist,l+1,Tbis,k)]] #niveau suivant
                else: #Ou alors on est arrivé au plus bas niveau et on a des feuilles
                    foret+=[[Teval,[]]]
        return foret


