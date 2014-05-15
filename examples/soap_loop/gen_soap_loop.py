from sps import *
import os
import profile
import sys
import pdb


logpath='/bell3/gqdong/statpot/results/mgl9_mgl10_mgl11_ml9_ml10_ml11/3xtop1_rmsd_mean_+10xNativeSelection+top3_rmsd_mean_+0.001xlbslope6-20-wn/runs/49331--2.36250,-35358.180_-2.030+0.105_-2.363+0.187_0.000_0.000/'
print logpath
bm=pickle.load(open(logpath+'cvmodel.pickle'))
bestpars=bm['allresult'][0]['repna']
bestpar=bm['originalresult'][0][0]['bestpar']
looplist=['mgl'+str(i) for i in range(4,21)]
#pdb.set_trace()
#bm['model']['bmtype']['dslist']=['ml12']
#bm['model']['scorers'][-2]['refs']=[bm['model']['scorers'][-1]]
#bm['model']['scorers']=bm['model']['scorers'][:-2]
so=scorer(model=bm['model'])

print so.assess(bm['allresult'][0]['bestpar'],slevel='top1_rmsd_mean_')
print so.assess_model()
print so.assess_ideal(slevel='top1_rmsd_mean_')
pdb.set_trace()
so.write_potential(filetype='hdf5')
so.calc_structurescore('/bell2/gqdong/rawdecoyfiles/tmp/')
print so.score[0]
pdb.set_trace()
