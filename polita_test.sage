attach("polITA.sage")

q1=Etat(1,1)
q2=Etat(2,2)
q3=Etat(3,2)

trans=[]
trans=trans+[[[(X1-1,1)],X1-1,q1,q2]]
trans=trans+[[[(X2-X1**2+1,1),(X1-1,-1)],X2-4,q2,q3]]
trans=trans+[[[(X2*X1-1,-1)],X2-1,q3,q1]]

ITA=PolITA([q1,q2,q3],[q1],[q3],trans)

