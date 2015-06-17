#(***************************************************************************************)
#(*                                                                                     *)
#(*                                                                                     *)
#(*                              HUOT Mathieu                                           *)
#(*                              GARNIER Remy                                           *)
#(*                    Licence 3 : stage de Mathématiques                               *)
#(*                           Version Q[X1,...Xn]                                       *)
#(*                             Parallelization                                         *)
#(*                                                                                     *)
#(*                                                                                     *)
#(***************************************************************************************)

#Active ou désactive la parallelization
def Parallelize_Mode(a):
    global ROOT_PAR
    global ROOT_PAR2
    global NORM_PAR
    global LIFT_PAR
    global PAR_PCC
    global PAR_PCC2
    ROOT_PAR=a
    ROOT_PAR2=a
    NORM_PAR=a
    LIFT_PAR=a
    PAR_PCC=a
    PAR_PCC2=a
    if a:
        print("parallelization switched on.")
    else:
        print("parallelization switched off.")
        
ROOT_PAR=True
@parallel
def RootPar(l,T,Pl,pl,P,p,i):
    return (RootCoding(l,T,Pl,pl,P,p),i)
    
ROOT_PAR2=True    
@parallel
def RootPar2(l,T,Pl,pl,P,p,i,j):
    return (RootCoding(l,T,Pl,pl,P,p),i,j)

NORM_PAR=True
@parallel
def NormalizePar(l,T,P,i):
    return (Normalize(l,T,P),i)

LIFT_PAR=True
@parallel
def LiftPar(PPtot,PPlist,l,T,k,i):
    return (Lift(PPtot,PPlist,l,T,k),i)

@parallel
def EvalPar(L,T,l,PPlist,i)
    return (Eval(L,T,l,PPlist)[0],i)

PAR_PCC=True
@parallel
def ParPreCalculCompleting(l,T,PP,i,j):
    return PreCalculCompleting(l,T,PP,i,j)
    
PAR_PCC2=True    
@parallel
def ParPreCalculCompleting2(l,T,PP,i):
    return PreCalculCompleting2(l,T,PP,i)
