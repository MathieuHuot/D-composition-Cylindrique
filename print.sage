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

def PrintArb2(Arb):
    def Aux2(Ar,l,s):
        if Ar==[]:
            return 0
        else:
            if l!=0:
                a=Espace(l-3)
                print(a+(''.join([str(_) for _ in s]))+" %s" %Ar[0])
            m=len(Ar[1])
            for j in range(m):
                Aux2(Ar[1][j],l+3,s+[j+1]+['.'])
            return 0
    Aux2(Arb,0,['-'])

def Aux(tree,level,number,max_depth):
        if tree==[]:
            return 0
        else:
            if level!=0:
                a=Espace(3*(level-1))
                print(a+(''.join([str(_) for _ in number]))+" %s" %tree[0])
            if level<max_depth:
                m=len(tree[1])
                for j in range(m):
                    Aux(tree[1][j],level+1,number+[j+1]+['.'],max_depth)
            return 0

def PrintArb(Arb,prof_max=2): #prof_max à changer)
    prompt = '> '
    print(" Bonjour\n Tapez 1 pour voir un lifting niveau par niveau\n et 2 pour avoir le lifting complet")
    n = int(raw_input(prompt))
    if n==2:
        PrintArb2(Arb)
    elif n!=1:
        print ("Vous n'avez pas suivi la consigne !")
    else:
        def PrintArbAux(Arb,i,s,count):
            print ("---%d niveau du lifting---" %i)
            print ("--------------------------")
            Aux(Arb,i,s,i+1)
            if count==0:
                print ("Entrez le numéro de la branche que vous souhaitez voir")
                num=""
                while num=="":
                    num = raw_input(prompt)
            else:
                print ("Entrez le numéro de la branche que vous souhaitez voir, ou back pour revenir en arrière")
                num=""
                while num=="":
                    num = raw_input(prompt)
                    if num=="back":
                        raise ValueError(i) 
            num=int(num)
            if num>len(Arb[1]):
                num=len(Arb[1])
            elif num<1:
                num=1
            if i<prof_max:
                try:
                    Arb2=Arb[1][num-1]
                    PrintArbAux(Arb2,i+1,s+[num]+['.'],1)
                except ValueError:
                    if i==1:
                        PrintArbAux(Arb,1,s,0)
                    else:
                        PrintArbAux(Arb,i-1,s,1)
        PrintArbAux(Arb,1,['-'],0)
