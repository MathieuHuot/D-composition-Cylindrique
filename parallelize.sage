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

PARALLELIZE=True

def Parallelize_Mode():
    if PARALLELIZE:
        PARALLELIZE=False
        print("parallelization switched off.")
    else:
        PARALLELIZE=True
        print("parallelization switched on.")

ROOT_PAR=True
@parallel
def RootPar(l,T,Pl,pl,P,p,i):
    return (RootCoding(l,T,Pl,pl,P,p),i)
    
ROOT_PAR2=True    
@parallel
def RootPar2(l,T,Pl,pl,P,p,i,j):
    return (RootCoding(l,T,Pl,pl,P,p),i,j)

ROOT_PAR3=True
@parallel
def NormalizePar(l,T,P,i):
    return (Normalize(l,T,P),i)

LIFT_PAR=True
@parallel
def LiftPar(PPtot,PPlist,l,T,k,i):
    return (Lift(PPtot,PPlist,l,T,k),i)

PAR_PCC=True
@parallel
def ParPreCalculCompleting(l,T,PP,i,j):
    return PreCalculCompleting(l,T,PP,i,j)
    
PAR_PCC2=True    
@parallel
def ParPreCalculCompleting2(l,T,PP,i):
    return PreCalculCompleting2(l,T,PP,i)
