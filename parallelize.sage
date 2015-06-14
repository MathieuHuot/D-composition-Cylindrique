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

@parallel
def RootPar(l,T,Pl,pl,P,p,i):
    return (RootCoding(l,T,Pl,pl,P,p),i)
    
@parallel
def RootPar2(l,T,Pl,pl,P,p,i,j):
    return (RootCoding(l,T,Pl,pl,P,p),i,j)

@parallel
def NormalizePar(l,T,P,i):
    return (Normalize(l,T,P),i)

@parallel
def LiftPar(PPtot,PPlist,l,T,k,i):
    return (Lift(PPtot,PPlist,l,T,k),i)

@parallel
def ParPreCalculCompleting(l,T,PP,i,j):
    return PreCalculCompleting(l,T,PP,i,j)
    
@parallel
def ParPreCalculCompleting2(l,T,PP,i):
    return PreCalculCompleting2(l,T,PP,i)
