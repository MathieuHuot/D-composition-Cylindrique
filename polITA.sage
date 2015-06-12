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

#Donne le nombre de variables d'un polynome
def nbVariables(P,maxi=4):
    if P in QQ:
        return 0
    while P in TdA[maxi]:
        maxi=maxi-1
    return maxi

#Fusion de deux ensembles en Ã©liminant les doublons
def Fusion(e1,e2):
    return list(set(e1+e2))

#Fait la liste des polynomes apparaissant dans un polITA
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

#Teste si un polynome vérifie une condition sur une cellule fixée
def Test(Con,cel):
    Cell=yolo[conv_lis_str(cel)]
    for po in Cell[1]:
        if po==Con:
            return True
    return False

#Décrit si un etat est accessible dans un ITA donné
def accessible(etat,ITA):
    Polist=listepol(ITA)
    print(Polist)
    EPolist=Elim(Polist)
    acc=[]
    qo=(ITA.initial)[0]
    l=qo.clock
    Access(EPolist,Polist,1,l,[1])
    print(yolo)
    etats=ITA.etats
    a=[1]
    Pere=yolo[conv_lis_str(a)]
    i=1
    while i <=l: 
        trouve=False
        for j in range(Pere[0]):
            Frere=yolo[conv_lis_str(a+[j])]
            for Co in Frere[1]:
                if Co[0]==TdV[i] and Co[1]==0:
                    a=a+[j]
                    trouve=True
                    Access(EPolist,Polist,i,i,a)
                    break
            if trouve:
                break
        Pere=Frere
        i=i+1
    acc=[Config(a,qo)]  #la config initiale 
    newacc=acc #Nouveaux états à parcourir
    oldacc=[]#Etats accessibles precedemment
    while acc != oldacc:
        for conf in newacc :
            if conf.etat==etat:
                return True
        ajout=[]
        for conf in newacc:
            if not (conf in oldacc):
                confAcc=Transition(conf,EPolist,Polist,ITA) # étatx accessibles en une étape
                ajout=ajout+confAcc
        oldacc=acc
        acc=Fusion(newacc,acc)
        newacc=ajout
    return False

#Donne la liste des configurations accessibles en une étape:
def Transition(conf,EPolist,Polist,ITA):
    q1=conf.etat
    cel=conf.cellule
    Tr=ITA.transitions
    etats=ITA.etats
    
    #On ajoute l'état suivant par passage du temps
    hauteur=len(list(cel))
    rang=cel[hauteur-1]
    ap=[cel[i] in range(hauteur-1)]
    pere=yolo[conv_lis_str(ap)]
    if rang<pere[0]:
        confAtteinte=[Config(ap+[rang+1],q2)]
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
                    newConf=Config(NewCel,q2)
                    confAtteinte=confAtteinte+[newConf]
    return confAtteinte

#Renvoie la cellule associée  après la transition q1->q2 en partant de cel
def AddCel(EPolist,Polist,cel,q1,q2,P,ITA):
    if q1.clock>q2.clock:
        #Cas où on remonte d'un cran:
        b=[cel[i] for i in range(q2.clock())]
        return b
        
    #Les autres cas:
    a=[cel[i] in range(len(a)-1)]
    Pere=yolo[conv_lis_str(a)]
    i=q1.clock()
    while i <=q2.clock(): 
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
        P=TdV[i+1]
        i=i+1
    return a
