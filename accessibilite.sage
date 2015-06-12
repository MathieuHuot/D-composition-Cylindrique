#(***************************************************************************************)
#(*                                                                                     *)
#(*                                                                                     *)
#(*                              HUOT Mathieu                                           *)
#(*                              GARNIER Remy                                           *)
#(*                    Licence 3 : stage de Mathématiques                               *)
#(*                           Version Q[X1,...Xn]                                       *)
#(*                            Lifting amélioré                                         *)
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



#Création du dictionnaire pour le lifting
yolo=dict()
yolo['1.']=[0,[],[]]

#Conversion de liste vers string
def conv_lis_str(lis):
    s= ''.join([str(_)+'.' for _ in lis])
    return s

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





#Fonction Lifting remplissant le dictionnaire
def Access(PolElim,PolIni,l,k,a): #Construction récursive de chaque niveau
    Cel=yolo[conv_lis_str(a)]
    if Cel[0]==0:

        Tsup=Cel[1]
        T=Cel[2]
        L,PP=LinePartition(PolElim[l-1],l,T)
        lon=len(PolIni[l-1])
        if L==[]: #Aucun polynome n'a de racine
            Tbis=T+[[1,TdV[l-1],1]] #X_l devient représentant de la ligne rÃ©elle
            Teval=[]
            for j in range(lon):
                P=PolIni[l-1][j]
                p=P.degree()
                s=Sign(l-1,T,P[p])  #Le signe d'un polynome sans racine réelle
                Teval=Teval+[(P,s)] #est celui de son coefficient dominant
            b=a+[0]
            yolo[conv_lis_str(b)]=[0,Tsup+Teval,Tbis]  
            if l<k:   #Si il reste un niveau à construire, on appelle récursivement
                Access(PolElim,PolIni,l+1,k,b)
            return ()
        else:
            foret=[]  #La ligne réelle est scindée par des racines de polynomes
            eval=[]   #On appelle donc completing pour avoir un représentant de
            L=Completing(l,T,L,PP) #de chaque cellule
            NewCel=[len(L),Tsup,T]
            yolo[conv_lis_str(a)]=NewCel
            for i in range(len(L)):
                Teval=[]
                ind=L[i][0] #L'indice i d'un Pi tel que L[i] code une racine de Pi 
                P=L[i][ind][2]
                r=L[i][ind][0]
                Tbis=T+[(r,P,Degree(l,T,P))] #Qu'on ajoute au système triangulaire
                for j in range(lon):
                    Pol=PolIni[l-1][j]
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
                b=a+[i]
                EvalP=Teval+Tsup
                yolo[conv_lis_str(b)]=[0,EvalP,Tbis]    
                if l<k: #On appelle récursivement sur chaque noeud la construction du 


                    Access(PolElim,PolIni,l+1,k,b) #niveau suivant
            return ()
    else:
        if l<k:
            for i in range(Cel[0]):
                Access(PolElim,PolIni,l+1,k,a+[i])
        return ()
