attach("Accessibilite")

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
    
class PolITA(Etat):
    """ Classe définissant un automate temporisé à garde polynomiale définie par:
    -Une liste d'états
    -Un dictionnaire de transitions """

    def __init__(self,etats,initials,finals,ttransitions):
        self.etats=etats
        self.initial=initials
        self.final=finals
        self.transitions=transitions

#Donne le nombre de variables d'un polynome
def nbVariables(P,max=10):
    while P in TdA[max]:
        max=max-1
    return max

#Fusion de deux ensembles
def Fusion(e1,e2):
    return list(set(e1+e2))

#Fait la liste des polynomes apparaissant dans un polITA
def listepol(ITA,max=10):
    L=[[TdV[i]] for i in range(max)]
    tr=ITA.transitions()
    for i in range(len(tr)):
        for P in tr[0]:
            j=nbVariables(P)
            L[j]=L[j]+[P]
    return L

#Teste si un polynome vérifie une condition sur une cellule fixée
def Test(Con,cel):
    Cell=yolo[conv_lis_str(cel)]
    for po in Cell[1]:
        if po=Con:
            return True
    return False

#Décrit si un etat est accessible dans un ITA donné
def access(etat,ITA):
    Polist=listepol(ITA)
    acc=ITA.initial()# Ajouter la cellule initiale 
    #TODO
    newacc=acc #Nouveaux états à parcourir
    oldacc=[]#Etats accessibles precedemment
    while acc != oldacc:
        for conf in newacc :
            if conf.etat==etat:
                return True
        ajout=[]
        for conf in newacc:
            if not (conf in oldacc):
                confAcc=Transition(conf,Polist,ITA) # etat accessibles en une étape
                ajout=ajout+confAcc
        oldacc=acc
        acc=Fusion(newAcc,acc)
        newAcc=ajout
    return False

#Donne la liste des configurations accessibles en une étape:
def Transition(conf,Polist,ITA):
    q1=conf.etat()
    cel=conf.cellule()
    Tr=ITA.transitions()
    etats=ITA.etats()
    etatAtteint=[]
    for q2 in etats:
        trans=Tr[(q1,q2)]
        if trans == 0:
            valide=True
            condition=trans[0]
            nbc=len(condition)
            i=0
            while valide and i<nbc:
                valide=Test(Condition[i],cel) #Teste si le polynome vérifie la condition 
                i=i+1
            if valide:
                NewCel=AddCel(cel,q1,q2,ITA)
    return




