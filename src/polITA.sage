#(***************************************************************************************)
#(*                                                                                     *)
#(*                                                                                     *)
#(*                              GARNIER Remy                                           *)
#(*                              HUOT Mathieu                                           *)
#(*                    Licence 3 : stage de MathÃ©matiques                              *)
#(*                            Version Q[X1,...Xn]                                      *)
#(*                           Accessibility decision                                    *)
#(*                                                                                     *)
#(*                                                                                     *)
#(***************************************************************************************)

import pdb
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
#OUTPUT: m integer : number of variables that appear in P
def nbVariables(P,maxi=4):
    if P in QQ:
        return 0
    while P in TdA[maxi]:
        maxi=maxi-1
    return maxi

#INPUT : e1 list
#        e2 list
#OUTPUT: e  list : fusion of e1 and e2 by removing repetitions
def Fusion(e1,e2):
    return list(set(e1+e2))

#INPUT : ITA PolITA
#OUTPUT: L   Q[X1,...,Xn] list : polynomials that appear in ITA
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
         if P!="None":
             j=nbVariables(P)
             L[j]+=[P]
    return L

#INPUT : Con int list
#        cel cell
#OUTPUT: b   boolean : tests if there exists a polynomial satisfying Con condition in cel
def Test(Con,cel,Arb):
    Cell=Arb[conv_lis_str(cel)]
    for po in Cell[1]:
        if po==Con:
            return True
    return False

#INPUT : etat state
#        ITA  PolITA
#OUTPUT: b boolean : True iff etat is accessible in the PolITA ITA
def accessible(etat,ITA):
    Arb=init_dict() #Cylindrical Decomposition Tree
    Polist=listepol(ITA)
    hmax=len(Polist)
    EPolist=Elim(Polist)
    acc=[]
    qo=(ITA.initial)[0]
    l=qo.clock
    Access(EPolist,Polist,1,l,[1],Arb)
    etats=ITA.etats
    a=[1]
    Pere=Arb[conv_lis_str(a)]
    i=1
    while i <=l:
        trouve=False
        for j in range(Pere[0]):
            Frere=Arb[conv_lis_str(a+[j])]
            for Co in Frere[1]:
                if Co[0]==TdV[i-1] and Co[1]==0:
                    a=a+[j]
                    trouve=True
                    if i<hmax:
                        Access(EPolist,Polist,i+1,i+1,a,Arb)
                    break
            if trouve:
                break
        Pere=Frere
        i=i+1
    acc=[Config(qo,a)]  #initial configuration
    newacc=acc #New states to go
    oldacc=[]#previously accessible states
    while set(acc) != set(oldacc):
        for conf in newacc :
            if conf.etat==etat:
                return True
        ajout=[]
        for conf in newacc:
            if not (conf in oldacc):
                confAcc=Transition(conf,EPolist,Polist,ITA,Arb,hmax) #reacheable states in 1 step
                ajout=ajout+confAcc
        oldacc=acc
        newacc=ajout
        acc=Fusion(newacc,acc)
    return False

#INPUT : conf: an initial configuration
#        EPolist: output of elim phase
#        Polist: initial polynomials
#        ITA: a polITA
#OUTPUT: confAtteinte : list of reachable configurations in one step
def Transition(conf,EPolist,Polist,ITA,Arb,hmax):
    q1=conf.etat
    cel=conf.cellule
    Tr=ITA.transitions
    etats=ITA.etats
    #We add the next state by time elapsing
    hauteur=len(cel)
    rang=cel[hauteur-1]
    ap=[cel[i] for i in range(hauteur-1)]
    pere=Arb[conv_lis_str(ap)]
    if rang<(pere[0]-1): #We check we are not in the end of the line
        confAtteinte=[Config(q1,ap+[rang+1])]
        print("Nouvelle cellule= %s" %(ap+[rang+1]))
        print(q1.nom)
    else:
        confAtteinte=[]
    #Then configurations obtained after tran
    for q2 in etats:
        for trans in Tr:
            if trans[2]==q1 and trans[3]==q2: 
                valide=True
                conditions=trans[0]
                nbc=len(conditions)
                i=0
                while valide and i<nbc:
                    valide=Test(conditions[i],cel,Arb) #Tests if the polynomial satisfies the condition
                    i=i+1
                if valide:
                    Update=trans[1]
                    NewCel=AddCel(EPolist,Polist,cel,q1,q2,Update,ITA,Arb,hmax) #Returns the cell after the transition
                    print("Nouvelle cellule= %s" %NewCel)
                    print(q2.nom)
                    newConf=Config(q2,NewCel)
                    confAtteinte=confAtteinte+[newConf]

    return confAtteinte

#Returns the linked cell after the transition q1->q2 from cell
#INPUT : EPolist 
#        Polist
#        cel
#        q1
#        q2
#        P
#        ITA
#OUTPUT: b    
def AddCel(EPolist,Polist,cel,q1,q2,P,ITA,Arb,hmax):
    if q1.clock>q2.clock:
        #Case when we go a level down:
        b=[cel[i] for i in range(q2.clock())]
        return b
        
    #Other cases:
    a=[cel[i] for i in range(len(cel)-1)]
    Pere=Arb[conv_lis_str(a)]
    for i in range(q1.clock,q2.clock+1): 
        trouve=False
        print(i)
        if P=="None": #When there is no update
            a=cel
            if i<hmax:
                Access(EPolist,Polist,i+1,i+1,a,Arb)
            P=TdV[i]
            Pere=Arb[conv_lis_str(a)]
            
        else:
            for j in range(Pere[0]):
                Frere=Arb[conv_lis_str(a+[j])]
                for Co in Frere[1]:
                    if Co[0]==P and Co[1]==0:
                        a=a+[j]
                        trouve=True
                        if i <hmax:
                            Access(EPolist,Polist,i+1,i+1,a,Arb)
                        break
                if trouve:
                    break
            Pere=Frere
            P=TdV[i]
    return a


