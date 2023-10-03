from __future__ import print_function
from sps import *
import profile
import sys


scaledsp1={'type':'scaledsp','pdbset':'Otherpot','features':'soap_pep','genmethod':'fpd','pm':'','refs':[],'ratio':[1.0]}
#firedock
#define the benchmark scheme
bmtype={'type':'dsscore','dslist':['fpdb','fpdu'],'criteria':'3xtop1_rmsdallif_mean_+2xtop3_rmsdallif_mean_','combine':'scoresum','bm':''}


#define the final model
model1={'scorers':[scaledsp1],'bmtype':bmtype,'searches':[],'sm':['powell'],'cvk':20,'cvm':'parallel',
    'clustermethod':{'clusteringperc':0.999,'clustercutoff':2.5,'scoreratio':[1,1,0,3,1],'pickwhichone':'median'}}

so=scorer(model=model1)
#ov=[0.00425609,  0.04942771,  0.38919625,  1.        ]
res=so.assess_basescore(slevel='top1_rmsd_mean_')
print(res)
print(so.assess_ideal())
pdb.set_trace()

ara=np.zeros(30)
tnr=range(0,30)
for p in tnr:
    tcp=int(10**(p/10.0))
    ara[p]=so.assess_model(slevel='top'+str(tcp)+'__nonempty_rmsd10ORirmsd4FIRST')

print(res)

print(so.assess_ideal())
pdb.set_trace()
#so.assess_sascore()
#so.assess_randomized_score(ov)
#sys.exit(0)
