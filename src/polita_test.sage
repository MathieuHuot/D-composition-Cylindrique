attach("polITA.sage")

#Test 1: The example of the presentation  
q1=Etat(1,1)
q2=Etat(2,2)
q3=Etat(3,2)

trans=[]
trans=trans+[[[(X1**2-X1-1,-1)],"None",q1,q2]]
trans=trans+[[[(X2**2*(2*X1-1)-1,1)],"None",q2,q3]]
trans=trans+[[[(X1**2+X2-5,-1)],"None",q3,q1]]

ITA=PolITA([q1,q2,q3],[q1],[q3],trans)



#Test 2:
p1=Etat(1,1)
p2=Etat(2,2)
p3=Etat(3,2)

trans=[]
trans=trans+[[[(X1**2-X1-1,-1)],X1,p1,p2]]
trans=trans+[[[(X2-X1**2+1,1),(X1-1,-1)],X2-4,p2,p3]]
trans=trans+[[[(X2*X1-1,-1)],X2-1,p3,p1]]

ITA_2=PolITA([p1,p2,p3],[p1],[p3],trans)
