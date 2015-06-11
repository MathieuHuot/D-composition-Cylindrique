
attach("accessibilite.sage")

class Etat:
    def __init__(self,nom,clock):
        self.nom=nom
        self.clock= clock

class Config:
    """Une configuration d'un automate est un couple constituÃ©:
    -d'un etat actif
    -d'une cellule actif"""
    def __init__(self,etat,cellule):
        self.cellule=cellule
        self.etat=etat
    
class PolITA:
    """ Classe dÃ©finissant un automate temporisÃ© Ã  garde polynomiale dÃ©finie par:
    -Une liste d'Ã©tats
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
    return L

#Teste si un polynome vÃ©rifie une condition sur une cellule fixÃ©e
def Test(Con,cel):
    Cell=yolo[conv_lis_str(cel)]
    for po in Cell[1]:
        if po==Con:
            return True
    return False

#DÃ©crit si un etat est accessible dans un ITA donnÃ©
def accessible(etat,ITA):
    Polist=listepol(ITA)
    print(Polist)
    EPolist=Elim(Polist)
    acc=[]
    qo=(ITA.initial)[0]
    l=qo.clock
    Access(EPolist,Polist,1,1,[1])
    etats=ITA.etats
    a=[1]
    Pere=yolo[conv_lis_str(a)]
    i=0
    while i <=l: 
        trouve=False
        for j in range(Pere[0]):
            Frere=yolo[conv_lis_str(a+[j])]
            for Co in Frere[1]:
                if Co[0]==TdV[i] and Co[1]==0:
                    a=a+[j]
                    trouve=True
                    Access(EPolist,Polist,i,1,a)
                    break
            if trouve:
                break
        Pere=Frere
        i=i+1
    acc=[Config(a,qo)]  #la config initiale 
    newacc=acc #Nouveaux Ã©tats Ã  parcourir
    oldacc=[]#Etats accessibles precedemment
    while acc != oldacc:
        for conf in newacc :
            if conf.etat==etat:
                return True
        ajout=[]
        for conf in newacc:
            if not (conf in oldacc):
                confAcc=Transition(conf,EPolist,Polist,ITA) # etat accessibles en une Ã©tape
                ajout=ajout+confAcc
        oldacc=acc
        acc=Fusion(newAcc,acc)
        newAcc=ajout
    return False

#Donne la liste des configurations accessibles en une Ã©tape:
def Transition(conf,EPolist,Polist,ITA):
    q1=conf.etat
    cel=conf.cellule
    Tr=ITA.transitions
    etats=ITA.etats
    confAtteinte=[] # Ajouter l'Ã©tat suivant sur la ligne, toujours atteint
    for q2 in etats:
        for trans in Tr:
            if trans[2]==q1 and trans[3]==q2: 
                valide=True
                condition=trans[0]
                nbc=len(condition)
                i=0
                while valide and i<nbc:
                    valide=Test(Condition[i],cel) #Teste si le polynome vÃ©rifie la condition 
                    i=i+1
                if valide:
                    P=trans[1]
                    NewCel=AddCel(EPolist,Polist,cel,q1,q2,P,ITA)
                    newConf=Config(NewCel,q2)
                    confAtteinte=confAtteinte+[newConf]
    return confAtteinte

def AddCel(EPolist,Polist,cel,q1,q2,P,ITA):
    if q1.clock>q2.clock:
        #Cas ou on remonte d'un cran:
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
                    Access(EPolist,Polist,i,1,a)
                    break
            if trouve:
                break
        Pere=Frere
        P=TdV[i+1]
        i=i+1
    return a
