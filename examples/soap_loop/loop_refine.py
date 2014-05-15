
from SOAP.scorer import *
import os
import profile
import sys
import pdb
import numpy as np

def looprefine_eval(cv):
    sprefiner=sprefine(self,dslist='ml12',bm='loop.fast.1000')
    return sprefiner.get_task(cv)
    

runenv.runpriority=-10
#mjloop=['mgl'+str(i) for i in range(6,17)]
#decoysets(dslist=mjloop).build()

ref1={'type':'sf','features':'d30','sftype':'spline4','par':[1,2,8,15],'parvalue':[1.0,1,1,1]}

scaledsp1={'type':'scaledsp','pdbset':'X_2.2A','features':'a158as158d30','genmethod':'ss1',
           'ExcludeLocal':False,'pm':'npend','refs':[ref1]}

#define the benchmark scheme
bmtype={'type':'sprefine','dslist':'ml12','criteria':'bestrmsd','bm':'loop.fast.10'}

#define the search positions
search1={'object':ref1,'key':'parvalue','pos':[0,1,2,3],'InitialGenerator':
    {'type':'values','values':[[0.00750637,  0.02547501,  0.38271168,  1.08959042]]}}

#define the final model
model1={'scorers':[ref1],'bmtype':bmtype,'searches':[search1],'sm':['fmin_powell']
    ,'cvk':5,'cvm':'parallel'}

so=scorer(model=model1)
#so.sprefiner.assess_local()
#pdb.set_trace()
so.assess_model()
#spl=sps(modellist=[model1])
#spl.initialize()
#spl.cv()
#spl.log()
#profile.run('so.write_potential(type=\'lib\')')
#1000 fast 4.9A
#1000 slow_broad 4.4A
#1000 default 4.5A
#1000 sgmd 4.48,4.71
