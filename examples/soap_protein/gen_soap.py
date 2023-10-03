from __future__ import print_function
from sps import *
import os
import profile
import sys
import pdb


logpath='/bell3/gqdong/statpot/results/DecoysRUS_baker1/top1_rmsd_mean_+10xNativeSelection+cc/runs/54251-nan,nan_nan+nan_nan+nan_3.276_0.000__nan+nan_nan+nan_0.000_0.000 /'

print(logpath)
bm=pickle.load(open(logpath+'cvmodel.pickle', 'rb'))
#bestpars=bm['allresult'][0]['repna']
bestpar=bm['originalresult'][0][0]['bestpar']
#looplist=['mgl'+str(i) for i in range(4,21)]

#bm['model']['bmtype']['dslist']=['mgl15']
#bm['model']['scorers'][-2]['refs']=[bm['model']['scorers'][-1]]
#bm['model']['scorers']=bm['model']['scorers'][:-2]
so=scorer(model=bm['model'])


print(so.assess(bm['allresult'][0]['bestpar'],slevel='top1_rmsd_mean_'))
print(so.assess_model())
print(so.assess_ideal(slevel='top1_rmsd_mean_'))
pdb.set_trace()
so.write_potential(filetype='hdf5')


so.calc_structurescore('/bell2/gqdong/rawdecoyfiles/tmp/')
print(so.score[0])
pdb.set_trace()
