#(***************************************************************************************)
#(*                                                                                     *)
#(*                                                                                     *)
#(*                              GARNIER Remy                                           *)
#(*                              HUOT Mathieu                                           *)
#(*                    Licence 3 : stage de Mathématiques                               *)
#(*                           Version Q[X1,...Xn]                                       *)
#(*                             Parallelization                                         *)
#(*                                                                                     *)
#(*                                                                                     *)
#(***************************************************************************************)

NB_CORES=sage.parallel.ncpus.ncpus()
USED_CORES=NB_CORES

#Switch parallelization mode on/off and ask for the number or cores to use
def Parallelize_Mode():
    print("Please enter 1 if you want to switch parallelization mode ON and 0 else.")
    prompt='>'
    a=int(raw_input(prompt))
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
        print("Parallelization switched on.")
        print("Please enter the number or cores you want to use for parallelization")
        global USED_CORES
        b=int(raw_input(prompt))
        if b<1:
            b=1
        elif b>NB_CORES:
            b=NB_CORES
        USED_CORES=b
    else:
        print("Parallelization switched off.")
        
ROOT_PAR=True
@parallel(ncpus=USED_CORES)
def RootPar(l,T,Pl,pl,P,p,i):
    return (RootCoding(l,T,Pl,pl,P,p),i)
    
ROOT_PAR2=True    
@parallel(ncpus=USED_CORES)
def RootPar2(l,T,Pl,pl,P,p,i,j):
    return (RootCoding(l,T,Pl,pl,P,p),i,j)

NORM_PAR=True
@parallel(ncpus=USED_CORES)
def NormalizePar(l,T,P,i):
    return (Normalize(l,T,P),i)

LIFT_PAR=True
@parallel(ncpus=USED_CORES)
def LiftPar(PPtot,PPlist,l,T,k,i):
    return (Lift(PPtot,PPlist,l,T,k),i)

@parallel(ncpus=USED_CORES)
def EvalPar(L,T,l,PPlist,i,j):
    return (Eval(L,T,l,PPlist,i)[0],j)

PAR_PCC=True
@parallel(ncpus=USED_CORES)
def ParPreCalculCompleting(l,T,PP,i,j):
    return PreCalculCompleting(l,T,PP,i,j)
    
PAR_PCC2=True    
@parallel(ncpus=USED_CORES)
def ParPreCalculCompleting2(l,T,PP,i):
    return PreCalculCompleting2(l,T,PP,i)

@parallel(ncpus=USED_CORES)
def ParCalcul(l,T,R,r,P,p,Uplet):
    R,r=IntRem(l,T,R,r,P,p)
    return PmVPol(l,T,P,p,R,r)
