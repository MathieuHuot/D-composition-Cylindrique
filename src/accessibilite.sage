#(***************************************************************************************)
#(*                                                                                     *)
#(*                                                                                     *)
#(*                              GARNIER Remy                                           *)
#(*                              HUOT Mathieu                                           *)
#(*                    Licence 3 : stage de MathÃ©matiques                               *)
#(*                           Version Q[X1,...Xn]                                       *)
#(*                          On-the-fly algorithm                                       *)
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

#Creation of the dictionnary for lifting
def init_dict():
    Arb=dict()
    Arb['1.']=[0,[],[]]
    return Arb

#INPUT : lis int list
#OUTPUT: s   string   : elements of lis separated by points
def conv_lis_str(lis):
    s= ''.join([str(_)+'.' for _ in lis])
    return s

#INPUT : P   Q[X1,...,Xl]
#        rac integer * (int * int list * Q[X1,...,Xl]) list
#OUTPUT: i   integer : position of P in rac if present and 0 else
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

#INPUT : PolElim Q[X1,...,Xn] list : polynomials of elim phase
#        PolIni  Q[X1,...,Xn] list : initial polynomials 
#        l       integer  : current level
#        k       integer  : maximum level
#        a       int list : codes the cell being currently treated
#        yolo    dictionnary : tree of cylindrical decomposition 
#OUTPUT: None
#Note  : function of lifting fulling yolo dictionnary
def Access(PolElim,PolIni,l,k,a,yolo): #recursive constrcution of every level
    Cel=yolo[conv_lis_str(a)]
    if Cel[0]==0:
        Tsup=Cel[1]
        T=Cel[2]
        L,PP=LinePartition(PolElim[l-1],l,T)
        lon=len(PolIni[l-1])
        if L==[]: #No polynomial has a root
            Tbis=T+[[1,TdV[l-1],1]] #X_l becomes the representant of the real line
            Teval=[]
            for j in range(lon):
                P=PolIni[l-1][j]
                p=P.degree()
                s=Sign(l-1,T,P[p])  #The sign of a polynomial with no roots
                Teval=Teval+[(P,s)] #is the one of its leading coefficient
            b=a+[0]
            yolo[conv_lis_str(b)]=[1,Tsup+Teval,Tbis]  
            if l<k:   #if there is a level left to build, we call access recursively
                Access(PolElim,PolIni,l+1,k,b,yolo)
            return ()
        else:
            foret=[]  #Real line is split by roots of polynomials
            eval=[]   #We call completing to obtain a representant of every cell
            L=Completing(l,T,L,PP) 
            NewCel=[len(L),Tsup,T]
            yolo[conv_lis_str(a)]=NewCel
            for i in range(len(L)):
                Teval=[]
                ind=L[i][0] #The index j of P_j so that L[i] codes a root of Pj 
                P=L[i][ind][2]
                r=L[i][ind][0]
                Tbis=T+[(r,P,Degree(l,T,P))] #We add to the triangular system
                for j in range(lon):
                    Pol=PolIni[l-1][j]
                    pos=RechP(Pol,L[i])
                    if pos>0:
                        Teval=Teval+[(L[i][pos][2],L[i][pos][1][0])]
                    else:
                        fini=False #as long as we may simplify
                        trouve=True #if we found a simplification : one more loop
                        sP=1 #Sign of Pol
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
                        Teval+=[(PolIni[l-1][j],sP)]
                b=a+[i]
                EvalP=Teval+Tsup
                yolo[conv_lis_str(b)]=[0,EvalP,Tbis]    
                if l<k: #We make a recursive call on every node to build next level
                    Access(PolElim,PolIni,l+1,k,b,yolo)
            return ()
    else:
        if l<k:
            for i in range(Cel[0]):
                Access(PolElim,PolIni,l+1,k,a+[i],yolo)
        return ()
