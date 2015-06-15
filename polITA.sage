# -*-coding:Latin-1 -*

#(***************************************************************************************)
#(*                                                                                     *)
#(*                                                                                     *)
#(*                              HUOT Mathieu                                           *)
#(*                              GARNIER Remy                                           *)
#(*                    Licence 3 : stage de MathÃ©matiques                              *)
#(*                            Version Q[X1,...Xn]                                      *)
#(*                       décision de l'accessibilité                                   *)
#(*                                                                                     *)
#(*                                                                                     *)
#(***************************************************************************************)

attach("accessibilite.sage")

class Etat:
    def __init__(self,nom,clock):
        self.nom=nom
        self.clock= clock

class Config:
    """Une configuration d'un automate est un couple constitué:
    -d'un etat actif
    -d'une cellule actif"""
    def __init__(self,etat,cellule):
        self.cellule=cellule
        self.etat=etat
    
class PolITA:
    """ Classe définissant un automate temporisé à garde polynomiale définie par:
    -Une liste d'états
    -Un dictionnaire de transitions """

    def __init__(self,etats,initials,finals,transitions):
        self.etats=etats
        self.initial=initials
        self.final=finals
        self.transitions=transitions

#INPUT : P Q[X1,...,Xm]
#        m integer (optionnal)
#OUTPUT: m integer : nombre de variables effectivement dans P
def nbVariables(P,maxi=4):
    if P in QQ:
        return 0
    while P in TdA[maxi]:
        maxi=maxi-1
    return maxi

#INPUT : e1 list
#        e2 list
#OUTPUT: e  list : fusion de e1 et e2 en enlevant les doublons
def Fusion(e1,e2):
    return list(set(e1+e2))

#INPUT : ITA PolITA
#OUTPUT: L   Q[X1,...,Xn] list : liste des polynômes apparaissant dans ITA
def listepol(ITA):
    L=[]
    Tr=ITA.transitions
    nbV=-1
    for i in range(len(Tr)):
        for P in Tr[i][0]:
            j=nbVariables(P[0])
            for k in range(nbV+1,j+1):
                L=L+[[TdV[k]]]
            nbV=max(j,nbV)
            L[j]=L[j]+[P[0]]

    for i in range(len(Tr)):
         P = Tr[i][1]
         j=nbVariables(P)
         L[j]+=[P]
    return L

#INPUT : Con int list
#        cel cellule
#OUTPUT: b   boolean : test si il existe un polynome vérifiant la condition Con sur la cellule cel
def Test(Con,cel):
    Cell=yolo[conv_lis_str(cel)]
    for po in Cell[1]:
        if po==Con:
            return True
    return False

#INPUT : etat etat
#        ITA  PolITA
#OUTPUT: b boolean : dit si l'état etat est accessible dans le PolITA ITA
def accessible(etat,ITA):
    Polist=listepol(ITA)
    EPolist=Elim(Polist)
    acc=[]
    qo=(ITA.initial)[0]
    l=qo.clock
    Access(EPolist,Polist,1,l,[1])
    etats=ITA.etats
    a=[1]
    Pere=yolo[conv_lis_str(a)]
    i=1
    while i <=l:
        trouve=False
        for j in range(Pere[0]):
            Frere=yolo[conv_lis_str(a+[j])]
            for Co in Frere[1]:
                if Co[0]==TdV[i-1] and Co[1]==0:
                    a=a+[j]
                    trouve=True
                    Access(EPolist,Polist,i,i,a)
                    break
            if trouve:
                break
        Pere=Frere
        i=i+1
    acc=[Config(qo,a)]  #la configuration initiale
    newacc=acc #Nouveaux états à parcourir
    oldacc=[]#Etats accessibles précedemment
    while set(acc) != set(oldacc):
        for conf in newacc :
            if conf.etat==etat:
                return True
        ajout=[]
        for conf in newacc:
            if not (conf in oldacc):
                confAcc=Transition(conf,EPolist,Polist,ITA) # états accessibles en une étape
                ajout=ajout+confAcc
        oldacc=acc
        newacc=ajout
        acc=Fusion(newacc,acc)
    return False

#INPUT : conf: Une configuration initiame
#        EPolist: L'ensemble des polynomes obtenus par Elim
#        Polist: L'ensemble des polynomes initials
#        ITA: un polITA
#OUTPUT: confAtteinte : liste des configurations accessibles en une étape
def Transition(conf,EPolist,Polist,ITA):
    q1=conf.etat
    cel=conf.cellule
    Tr=ITA.transitions
    etats=ITA.etats
    print(cel)
    #On ajoute l'état suivant par passage du temps
    hauteur=len(cel)
    rang=cel[hauteur-1]
    ap=[cel[i] for i in range(hauteur-1)]
    pere=yolo[conv_lis_str(ap)]
    if rang<(pere[0]-1): #On vérifie que l'on est pas en bout de ligne
        confAtteinte=[Config(q1,ap+[rang+1])]
    else:
        confAtteinte=[]
    #Puis les configurations obtenues après tran
    for q2 in etats:
        for trans in Tr:
            if trans[2]==q1 and trans[3]==q2: 
                valide=True
                conditions=trans[0]
                nbc=len(conditions)
                i=0
                while valide and i<nbc:
                    valide=Test(conditions[i],cel) #Teste si le polynome vérifie la condition 
                    i=i+1
                if valide:
                    Update=trans[1]
                    NewCel=AddCel(EPolist,Polist,cel,q1,q2,Update,ITA) #Renvoie la cellule en fin de transition
                    newConf=Config(q2,NewCel)
                    confAtteinte=confAtteinte+[newConf]

    return confAtteinte

#Renvoie la cellule associée  après la transition q1->q2 en partant de cel
#INPUT : EPolist 
#        Polist
#        cel
#        q1
#        q2
#        P
#        ITA
#OUTPUT: b    
def AddCel(EPolist,Polist,cel,q1,q2,P,ITA):
    if q1.clock>q2.clock:
        #Cas où on remonte d'un cran:
        b=[cel[i] for i in range(q2.clock())]
        return b
        
    #Les autres cas:
    a=[cel[i] for i in range(len(cel)-1)]
    Pere=yolo[conv_lis_str(a)]
    i=q1.clock
    while i <=q2.clock: 
        trouve=False
        for j in range(Pere[0]):
            Frere=yolo[conv_lis_str(a+[j])]
            for Co in Frere[1]:
                if Co[0]==P and Co[1]==0:
                    a=a+[j]
                    trouve=True
                    Access(EPolist,Polist,i,i,a)
                    break
            if trouve:
                break
        Pere=Frere
        P=TdV[i]
        i=i+1
    return a

