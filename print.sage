#(***************************************************************************************)
#(*                                                                                     *)
#(*                                                                                     *)
#(*                              HUOT Mathieu                                           *)
#(*                              GARNIER Remy                                           *)
#(*                    Licence 3 : stage de Mathématiques                               *)
#(*                           Version Q[X1,...Xn]                                       *)
#(*                            Print lifting tree                                       *)
#(*                                                                                     *)
#(*                                                                                     *)
#(***************************************************************************************)

#INPUT : l integer
#OUTPUT: a string : un blanc de longueur l
def Espace(l):
    a=" "
    for i in range(l):
        a=a+"  "
    return a

#Renvoie une version indentée et numérotée de l'arbre construit par le Lifting
def PrintArb(Arb):
    def Aux(Ar,l,s):
        if Ar==[]:
            return 0
        else:
            if l!=0:
                a=Espace(l-3)
                print(a+(''.join([str(_) for _ in s]))+" %s" %Ar[0])
            m=len(Ar[1])
            for j in range(m):
                Aux(Ar[1][j],l+3,s+[j+1]+['.'])
            return 0
    Aux(Arb,0,['-'])
