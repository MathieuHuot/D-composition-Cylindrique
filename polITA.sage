
#attach("Accessibilite")



class Etat:
    def __init__(self,nom,clock):
        self.nom=nom
        self.clock= clock

class PolITA(Etat):
    """ Classe définissant un automate temporisé à garde polynomiale définie par:
    -Une liste d'états
    -Un dictionnaire de transitions """


    def __init__(self,etats,initials,finals,transitions):
        self.etats=etats
        self.initial=initials
        self.final=finals
        self.transitions=transitions



def nbVariables(P,max=10):
    while P in TdA[max]:
        max=max-1
    return max


