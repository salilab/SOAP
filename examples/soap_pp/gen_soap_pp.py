from sps import *
import os
import profile
import sys
import pdb


logpath='/bell3/gqdong/statpot/results/ppd4s/top1000_rlrank__rmsd10ORirmsd4FIRST+top1000_rlrank__rmsd5ORirmsd2FIRST/runs/39882-2.62460,-261.084_2.808+0.035_2.625+0.073_2.773_0.000/'
bm=pickle.load(open(logpath+'cvmodel.pickle'))
bestpars=bm['allresult'][0]['bestpar']#bm['allresult'][0]['repna']
bm['model']['bmtype']['dslist']=['ppd4']
bm['model']['scorers'][-2]['refs']=[bm['model']['scorers'][-1]]
bm['model']['scorers']=bm['model']['scorers']
#del bm['model']['bmtype']['filters']
bm['model']['searches']=[]
#pdb.set_trace()

so=scorer(model=bm['model'])
#bm['model']['scorers'][1:2]+

#print so.assess(bm['allresult'][0]['bestpar'])
print so.assess_model(slevel='top10__nonempty_rmsd10ORirmsd4FIRST')
#pdb.set_trace()
so.write_potential(filetype='hdf5')
so.calc_structurescore('/bell2/gqdong/rawdecoyfiles/ppd/temp/')
#print so.score[0]
pdb.set_trace()
