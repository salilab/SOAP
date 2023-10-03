from __future__ import print_function
from sps import *
import profile
import sys


scaledsp1={'type':'scaledsp','pdbset':'Otherpot','features':'soap_protein_od','genmethod':'','pm':'','refs':[],'ratio':[1.0]}

#define the benchmark scheme
looplist=['mgl'+str(i) for i in range(6,21)]
mjloop=['ml'+str(i) for i in range(4,13)]
bmtype={'type':'dsscore','dslist':['DecoysRUS','baker1'],'criteria':'top1_rmsd_mean_+10xNativeSelection+cc','combine':'scoresum','bm':''}
#+2xtop3_rmsd_mean_ top1_rmsd_mean_+10x +cc
#define the final model ,'casp58','baker1'
model1={'scorers':[scaledsp1],'bmtype':bmtype,'searches':[],'sm':['powell'],'cvk':20,'cvm':'parallel',
    'clustermethod':{'clusteringperc':0.999,'clustercutoff':2.5,'scoreratio':[1,1,0,3,1],'pickwhichone':'median'}}

so=scorer(model=model1)
#ov=[0.00425609,  0.04942771,  0.38919625,  1.        ]
res=so.assess_basescore()
print(res)
print(so.assess_ideal())
pdb.set_trace()
#so.assess_sascore()
#so.assess_randomized_score(ov)
#sys.exit(0)
